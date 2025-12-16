use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};

use anyhow::{anyhow, Context, Result};
use bio::alphabets::dna::revcomp;
use block_aligner::scan_block::{Block, PaddedBytes};
use block_aligner::scores::{Gaps, NucMatrix};
use clap::Parser;
use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use seq_io::fastq::{Reader, Record};
use serde::Serialize;

#[derive(Parser, Debug)]
#[command(author, version, about = "Multiplexed amplicon demultiplexing and quantification")]
struct Args {
    #[arg(long)]
    r1: PathBuf,
    #[arg(long)]
    r2: PathBuf,
    #[arg(long)]
    samples: PathBuf,
    #[arg(long)]
    genes: PathBuf,
    #[arg(long)]
    output: PathBuf,
    #[arg(long, default_value_t = num_cpus::get())]
    threads: usize,
    #[arg(long, default_value_t = false)]
    debug: bool,
}

#[derive(Debug, Clone)]
struct SampleRecord {
    id: String,
    primer_a_core: String,
    primer_b_core: String,
    primer_a_core_rc: String,
    primer_b_core_rc: String,
}

#[derive(Debug, Clone)]
struct GeneRecord {
    id: String,
    wt: String,
}

#[derive(Serialize)]
struct OutputRow {
    sample_id: String,
    gene_id: String,
    r#type: String,
    sequence: String,
    count: u64,
    percentage: f64,
}

#[derive(Copy, Clone, Debug)]
enum Orientation {
    Forward,
    Swapped,
}

fn read_genes(path: &PathBuf) -> Result<Vec<GeneRecord>> {
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(true)
        .flexible(true)
        .from_path(path)?;
    let mut genes = Vec::new();
    for rec in rdr.records() {
        let rec = rec?;
        let id = rec
            .get(0)
            .ok_or_else(|| anyhow!("Gene_ID missing"))?
            .to_string();
        let seq = rec
            .get(1)
            .ok_or_else(|| anyhow!("WT_Sequence missing"))?
            .to_ascii_uppercase();
        genes.push(GeneRecord { id, wt: seq });
    }
    Ok(genes)
}

// gene_seeds: 앞부분 k=10bp를 seeds로 사용, 프라이머에 포함된 오버행을 자동으로 제거
fn strip_overhang_auto(primer: &str, gene_seeds: &[String]) -> String {
    let p = primer.to_ascii_uppercase();
    let mut best = p.len();
    for seed in gene_seeds {
        if let Some(pos) = p.find(seed) {
            if pos < best {
                best = pos;
            }
        }
    }
    if best < p.len() {
        p[best..].to_string()
    } else {
        p
    }
}

fn read_samples(path: &PathBuf, gene_seeds: &[String]) -> Result<Vec<SampleRecord>> {
    let mut rdr = csv::Reader::from_path(path)?;
    let mut out = Vec::new();
    for rec in rdr.records() {
        let rec = rec?;
        let sample_id = rec
            .get(0)
            .ok_or_else(|| anyhow!("Sample_ID missing"))?
            .to_string();
        let primer_a_raw = rec
            .get(1)
            .ok_or_else(|| anyhow!("Primer_A missing"))?
            .to_string();
        let primer_b_raw = rec
            .get(2)
            .ok_or_else(|| anyhow!("Primer_B missing"))?
            .to_string();

        let primer_a_core = strip_overhang_auto(&primer_a_raw, gene_seeds);
        let primer_b_core = strip_overhang_auto(&primer_b_raw, gene_seeds);

        out.push(SampleRecord {
            id: sample_id,
            primer_a_core_rc: String::from_utf8(revcomp(primer_a_core.as_bytes()))?,
            primer_b_core_rc: String::from_utf8(revcomp(primer_b_core.as_bytes()))?,
            primer_a_core,
            primer_b_core,
        });
    }
    Ok(out)
}

fn load_fastq_pairs(r1: &PathBuf, r2: &PathBuf) -> Result<Vec<(Vec<u8>, Vec<u8>)>> {
    let r1_file = File::open(r1).with_context(|| format!("Opening {:?}", r1))?;
    let r2_file = File::open(r2).with_context(|| format!("Opening {:?}", r2))?;
    let mut r1_reader = Reader::new(MultiGzDecoder::new(BufReader::new(r1_file)));
    let mut r2_reader = Reader::new(MultiGzDecoder::new(BufReader::new(r2_file)));

    let mut pairs = Vec::new();
    loop {
        let r1_rec = match r1_reader.next() {
            Some(res) => res?,
            None => break,
        };
        let r2_rec = match r2_reader.next() {
            Some(res) => res?,
            None => break,
        };
        pairs.push((r1_rec.seq().to_vec(), r2_rec.seq().to_vec()));
    }
    Ok(pairs)
}

// Fuzzy primer search
fn find_primer_offset_fuzzy(read: &[u8], primer: &[u8], window: usize, max_mismatch: usize) -> Option<usize> {
    let end = read.len().min(window.saturating_add(primer.len()));
    for start in 0..=end.saturating_sub(primer.len()) {
        let mism = read[start..start + primer.len()]
            .iter()
            .zip(primer.iter())
            .filter(|(a, b)| a != b)
            .count();
        if mism <= max_mismatch {
            return Some(start);
        }
    }
    None
}

fn detect_sample<'a>(
    samples: &'a [SampleRecord],
    r1: &[u8],
    r2: &[u8],
) -> Option<(&'a SampleRecord, Orientation)> {
    const WIN: usize = 200;
    const MM: usize = 3;
    for s in samples {
        let a = s.primer_a_core.as_bytes();
        let b = s.primer_b_core.as_bytes();
        let a_rc = s.primer_a_core_rc.as_bytes();
        let b_rc = s.primer_b_core_rc.as_bytes();

        if find_primer_offset_fuzzy(r1, a, WIN, MM).is_some()
            && find_primer_offset_fuzzy(r2, b, WIN, MM).is_some()
        {
            return Some((s, Orientation::Forward));
        }
        if (find_primer_offset_fuzzy(r1, b, WIN, MM).is_some()
            || find_primer_offset_fuzzy(r1, b_rc, WIN, MM).is_some())
            && (find_primer_offset_fuzzy(r2, a, WIN, MM).is_some()
                || find_primer_offset_fuzzy(r2, a_rc, WIN, MM).is_some())
        {
            return Some((s, Orientation::Swapped));
        }
    }
    None
}

/// Strip primers from both ends (core primers)
fn trim_primers(mut seq: Vec<u8>, sample: &SampleRecord) -> Vec<u8> {
    let pa = sample.primer_a_core.as_bytes();
    let pb = sample.primer_b_core.as_bytes();
    let pa_rc = sample.primer_a_core_rc.as_bytes();
    let pb_rc = sample.primer_b_core_rc.as_bytes();

    if seq.starts_with(pa) {
        seq = seq[pa.len()..].to_vec();
    } else if seq.starts_with(pa_rc) {
        seq = seq[pa_rc.len()..].to_vec();
    }
    if seq.ends_with(pb) {
        seq.truncate(seq.len() - pb.len());
    } else if seq.ends_with(pb_rc) {
        seq.truncate(seq.len() - pb_rc.len());
    }
    seq
}

/// Global alignment via block-aligner (완화된 기준)
fn best_alignment<'a>(genes: &'a [GeneRecord], seq: &[u8]) -> (&'a str, bool, f64) {
    if genes.is_empty() {
        return ("Unknown", false, f64::NEG_INFINITY);
    }
    let max_ref = genes.iter().map(|g| g.wt.len()).max().unwrap_or(0).max(1);
    let max_query = seq.len().max(1);
    let matrix = NucMatrix::new_simple(2, -2);
    let gaps = Gaps { open: -4, extend: -1 };
    let max_len = max_ref.max(max_query);
    let max_block_size = max_len.next_power_of_two().max(32).min(16_384);
    let min_block_size = 32;
    let mut block = Block::<true, false>::new(max_query, max_ref, max_block_size);

    let q = PaddedBytes::from_bytes::<NucMatrix>(seq, max_block_size);

    let mut best_id = "Unknown";
    let mut best_norm = f64::NEG_INFINITY;
    let mut second_norm = f64::NEG_INFINITY;
    let mut is_wt = false;

    for g in genes {
        let r = PaddedBytes::from_bytes::<NucMatrix>(g.wt.as_bytes(), max_block_size);
        block.align(&q, &r, &matrix, gaps, min_block_size..=max_block_size, 50);
        let res = block.res();
        let score = res.score as f64;
        let denom = 2.0 * (seq.len().max(g.wt.len())) as f64;
        let norm = score / denom;
        if norm > best_norm {
            second_norm = best_norm;
            best_norm = norm;
            best_id = g.id.as_str();
            let length_close = (seq.len() as isize - g.wt.len() as isize).abs() <= 10;
            let near_perfect = norm >= 0.90;
            is_wt = length_close && near_perfect;
        } else if norm > second_norm {
            second_norm = norm;
        }
    }

    if (best_norm - second_norm) < 0.001 {
        ("Unknown", false, best_norm)
    } else {
        (best_id, is_wt, best_norm)
    }
}

// Merge; fallback to concatenation if overlap fails.
fn merge_reads(r1: &[u8], r2: &[u8], min_overlap: usize, max_mismatch: usize) -> Option<Vec<u8>> {
    let max_olap = usize::min(r1.len(), r2.len());
    for olap in (min_overlap..=max_olap).rev() {
        let mismatches = (0..olap)
            .filter(|&i| r1[r1.len() - olap + i] != r2[i])
            .count();
        if mismatches <= max_mismatch {
            let mut merged = Vec::with_capacity(r1.len() + r2.len() - olap);
            merged.extend_from_slice(r1);
            merged.extend_from_slice(&r2[olap..]);
            return Some(merged);
        }
    }
    None
}

fn main() -> Result<()> {
    let args = Args::parse();
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .context("Configuring rayon")?;

    // 먼저 gene을 읽고 seeds 계산
    let genes = read_genes(&args.genes).context("Reading genes")?;
    let max_ref_len = genes.iter().map(|g| g.wt.len()).max().unwrap_or(0);
    let gene_seeds: Vec<String> = genes
        .iter()
        .map(|g| g.wt.chars().take(10).collect::<String>())
        .collect();
    let samples = read_samples(&args.samples, &gene_seeds).context("Reading samples")?;
    let gene_map: HashMap<String, GeneRecord> =
        genes.iter().cloned().map(|g| (g.id.clone(), g)).collect();
    let pairs = load_fastq_pairs(&args.r1, &args.r2).context("Loading FASTQ pairs")?;

    let total_pairs = AtomicUsize::new(0);
    let demux_hits = AtomicUsize::new(0);
    let merged_fallback = AtomicUsize::new(0);
    let primer_miss = AtomicUsize::new(0);

    // counts: key = (Sample, Gene, Type), value = (representative sequence, count)
    let counts: HashMap<(String, String, String), (String, u64)> = pairs
        .par_iter()
        .filter_map(|(r1_raw, r2_raw)| {
            let r1: Vec<u8> = r1_raw.iter().map(|b| b.to_ascii_uppercase()).collect();
            let r2: Vec<u8> = r2_raw.iter().map(|b| b.to_ascii_uppercase()).collect();

            total_pairs.fetch_add(1, Ordering::Relaxed);
            let (sample, orientation) = match detect_sample(&samples, &r1, &r2) {
                Some(x) => x,
                None => {
                    primer_miss.fetch_add(1, Ordering::Relaxed);
                    return None;
                }
            };
            demux_hits.fetch_add(1, Ordering::Relaxed);

            let merged = match merge_reads(&r1, &r2, 15, 4) {
                Some(m) => m,
                None => {
                    merged_fallback.fetch_add(1, Ordering::Relaxed);
                    let mut v = Vec::with_capacity(r1.len() + r2.len());
                    v.extend_from_slice(&r1);
                    v.extend_from_slice(&r2);
                    v
                }
            };
            let oriented = match orientation {
                Orientation::Forward => merged,
                Orientation::Swapped => revcomp(&merged),
            };

            // 길이 과다 방지: 참조 길이 + 30으로 캡
            let mut oriented = oriented;
            let max_len_cap = max_ref_len + 30;
            if oriented.len() > max_len_cap {
                oriented.truncate(max_len_cap);
            }

            let trimmed = trim_primers(oriented, sample);
            let (gene_id, is_wt, _norm) = best_alignment(&genes, &trimmed);

            // Representative sequence for aggregation
            let seq_rep = if is_wt {
                gene_map
                    .get(gene_id)
                    .map(|g| g.wt.clone())
                    .unwrap_or_else(|| String::from_utf8_lossy(&trimmed).into_owned())
            } else {
                let g_len = gene_map
                    .get(gene_id)
                    .map(|g| g.wt.len())
                    .unwrap_or(trimmed.len());
                let mut s = trimmed.clone();
                if s.len() > g_len {
                    s.truncate(g_len);
                }
                String::from_utf8_lossy(&s).into_owned()
            };

            let variant_type = if is_wt { "WT" } else { "Mutant" }.to_string();
            Some(((sample.id.clone(), gene_id.to_string(), variant_type), seq_rep, 1u64))
        })
        .fold(HashMap::new, |mut acc, (key, seq_rep, c)| {
            let entry = acc.entry(key).or_insert((seq_rep, 0u64));
            entry.1 += c;
            acc
        })
        .reduce(HashMap::new, |mut a, b| {
            for (k, v) in b {
                let entry = a.entry(k).or_insert(v.clone());
                entry.1 += v.1;
                if entry.0.is_empty() {
                    entry.0 = v.0;
                }
            }
            a
        });

    // Totals per Sample_ID + Gene_ID (WT+Mutant 합)
    let mut totals: HashMap<(String, String), u64> = HashMap::new();
    for ((sample, gene, _), (_, c)) in counts.iter() {
        *totals.entry((sample.clone(), gene.clone())).or_insert(0) += *c;
    }

    let mut wtr = csv::Writer::from_path(&args.output)?;
    for ((sample, gene, vtype), (seq, count)) in counts.iter() {
        let total = *totals.get(&(sample.clone(), gene.clone())).unwrap_or(&0);
        let denom = if total == 0 { 1 } else { total };
        // 퍼센트는 샘플-유전자 내 WT/Mutant 비율
        let pct = (*count as f64) * 100.0 / (denom as f64);
        wtr.serialize(OutputRow {
            sample_id: sample.clone(),
            gene_id: gene.clone(),
            r#type: vtype.clone(),
            sequence: seq.clone(),
            count: *count,
            percentage: pct,
        })?;
    }
    wtr.flush()?;

    if args.debug {
        eprintln!(
            "DEBUG: total_pairs={}, demux_hits={}, primer_miss={}, merged_fallback={}",
            total_pairs.load(Ordering::Relaxed),
            demux_hits.load(Ordering::Relaxed),
            primer_miss.load(Ordering::Relaxed),
            merged_fallback.load(Ordering::Relaxed)
        );
        let mut gene_counts: HashMap<&str, u64> = HashMap::new();
        for ((_, g, _), (_, c)) in counts.iter() {
            *gene_counts.entry(g.as_str()).or_insert(0) += *c;
        }
        let mut gc: Vec<_> = gene_counts.into_iter().collect();
        gc.sort_by_key(|(_, v)| std::cmp::Reverse(*v));
        eprintln!(
            "DEBUG: gene counts (top 10): {:?}",
            gc.into_iter().take(10).collect::<Vec<_>>()
        );
    }

    Ok(())
}
