use bio::pattern_matching::myers::{long, Myers};
use std::str::from_utf8;
use seq_io::fasta::{Reader, Record};
use anyhow::{Result, bail};
use clap::Parser;
use std::collections::HashMap;
use regex::Regex;
use serde::Serialize;
use simple_logger::SimpleLogger;
use time::macros::format_description;

#[derive(Debug)]
enum MyersEnum {
    Long(Box<long::Myers<u8>>),
    Short(Box<Myers>),
}

#[allow(dead_code)]
#[derive(Debug)]
struct LeaderSequence {
    id: String,
    seq: Vec<u8>,
    len: usize,
    myers: MyersEnum,
}

// fh_out.write("sequence_id\tsequence\tsequence_alignment\tgermline_alignment\tfwr1\tcdr1\tfwr2\tcdr2\tfwr3\tcdr3\tfwr4\n")
#[derive(Debug, Serialize)]
struct OutRecord<'a> {
    sequence_id: &'a str,
    sequence: &'a str,
    sequence_alignment: &'a str,
    germline_alignment: &'a str,
    fwr1: &'a str,
    cdr1: &'a str,
    fwr2: &'a str,
    cdr2: &'a str,
    fwr3: &'a str,
    cdr3: &'a str,
    fwr4: &'a str,
}
fn load_leading(f: &str) -> Result<Vec<LeaderSequence>> {
    let (reader, _format) = match niffler::send::from_path(f) {
        Ok(x) => x,
        Err(e) => bail!("Error opening file for reading: {}: {}",f , e)
    };
    let mut reader = Reader::new(reader);

    let mut v = vec![];
    while let Some(record) = reader.next() {
        let record = record?;
 
        let id = from_utf8(record.head())?.to_string();
        // let seq = record.seq();
        let seq = record.seq_lines().fold(
                vec![],
                |mut acc, x| {
                    acc.extend(x.to_ascii_uppercase());
                    acc
                }
            );
        let len = seq.len();
        let myers = if len > 64 {
            MyersEnum::Long(Box::new(long::Myers::<u8>::new(seq.clone())))
        } else {
            MyersEnum::Short(Box::new(Myers::<u64>::new(seq.clone())))
        };
        v.push(
            LeaderSequence {
                id,
                seq,
                len,
                myers
            }
        )
    }

    Ok(v)
}

fn load_annot(f: &str) -> Result<HashMap<String, usize>> {
    let re = Regex::new(r"\(\d+\):\((\d+)\-(\d+)\):\(\d+\-\d+\)").unwrap();
    let (reader, _format) = match niffler::send::from_path(f) {
        Ok(x) => x,
        Err(e) => bail!("Error opening file for reading: {}: {}",f , e)
    };

    let mut reader = Reader::new(reader);
    let mut annot = HashMap::new();
    for record in reader.records() {
        let record = record?;
        let id = record.id()?.to_string();
        let desc = match record.desc(){
            Some(d) => d?.to_string(),
            _ => "".to_string()
        };
        // >CGGACACAGTGCTGCC_1 801 29100.86 IGKV1-5*04(287):(241-506):(0-265):92.86,IGKV1-12*01(287):(241-530):(0-280):85.11,IGKV1D-12*02(287):(241-530):(0-280):85.11 * IGKJ2*03(39):(521-559):(0-38):97.44 IGKC*04(321):(559-800):(0-241):100.00 CDR1(319-336):88.89=CAGAGTATTAGTAATTGG CDR2(388-396):100.00=AAGGCATCT CDR3(502-531):100.00=TGCCAGCAATATAGTGACTTGTACAGTTTT
        let mut desc = desc.split_whitespace();
        let c_region = desc.nth(5);
        if let Some(c) = c_region {
            if c == "*"  {
                continue;
            }else{
                let start = re.captures_iter(c).map(|cap|  cap[1].parse::<usize>().unwrap()).next();
                if let Some(s)= start {
                    annot.insert(id, s);
                }
            }
        }
    }
    Ok(annot)
}

fn find_leading(myers: &MyersEnum, sequence: &[u8], max_dist: u8) -> Result<(usize, usize, usize)> {
    let occ: (usize, usize, usize) = match myers {
        MyersEnum::Long(m)=> {
            let mut myers = m.clone();
            let mut matches = myers.find_all_lazy(sequence, max_dist as usize);
            if let Some((best_end, dist)) = matches.by_ref().min_by_key(|&(_, dist)| dist) {
                if let Some((best_start, _)) = matches.hit_at(best_end){
                    (best_start, best_end + 1, dist)
                } else {
                    (0, 0, 0)
                }
                
            } else {
                (0, 0, 0)
            }
        },

        MyersEnum::Short(m)=> {
            let mut myers = m.clone();
            let mut matches = myers.find_all_lazy(sequence, max_dist);
            if let Some((best_end, dist)) = matches.by_ref().min_by_key(|&(_, dist)| dist) {
                if let Some((best_start, _)) = matches.hit_at(best_end){
                    (best_start, best_end + 1, dist as usize)
                } else {
                    (0, 0, 0)
                }
            } else {
                (0, 0, 0)
            }
        }
    };
    Ok(occ)
}

#[derive(Parser, Debug)]
#[command(version, about, long_about = None, arg_required_else_help(true))]
struct Args {
    /// leader fasta file
    #[arg(long, required(true),)]
    leader_ref: String,

    /// airr tsv file
    #[arg(long, required(true),)]
    airr: String,

    /// trust4 annot file
    #[arg(long, required(true),)]
    annot: String,

    /// max hamming distance
    #[arg(long, required(false), default_value_t=3)]
    dist: u8,

    /// output file
    #[arg(long, required(true),)]
    outfile: String,
}

fn main() -> Result<()> {
    let args = Args::parse();

    SimpleLogger::new()
        .with_local_timestamps()
        .with_timestamp_format(format_description!("[year]-[month]-[day] [hour]:[minute]:[second]"))
        .init()
        .unwrap();
    log::info!("Starting...");

    let leader_data = load_leading(&args.leader_ref)?;
    let c_data = load_annot(&args.annot)?;

    let (reader, _) = match niffler::send::from_path(&args.airr){
        Ok(x) => x,
        Err(e) => bail!("Error opening file for reading: {}: {}", args.airr , e)
    };
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(reader);
    let header = rdr.headers()?;
    let header_fields: HashMap<String, usize> = header.iter().enumerate().map(|(i, x)| (x.to_string(), i)).collect();
    let mut total = 0;
    let mut found = 0;
    let mut _perfect = 0;
    let mut wtr = match csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(&args.outfile) {
            Ok(wtr) => wtr,
            Err(e) => panic!("Error opening {} for writing: {}", &args.outfile, e),
        };
        // .from_writer(io::stdout());
    while let Some(Ok(record)) = &rdr.records().next() {
        total += 1;
        let sequence = record.get(header_fields["sequence"]).unwrap().as_bytes();
        let contig_name = record.get(header_fields["sequence_id"]).unwrap();
        let v_call = record.get(header_fields["v_call"]).unwrap();
        let pattern = format!(r"{}\b", regex::escape(v_call));
        let v_re = Regex::new(&pattern).unwrap();
        let cdr1 = record.get(header_fields["cdr1"]).unwrap().as_bytes();
        let leading_max = if !cdr1.is_empty() {
            sequence.windows(cdr1.len()).position(|window| window == cdr1).unwrap()
        } else {
            sequence.len()
        };
        let cdr2 = record.get(header_fields["cdr2"]).unwrap().as_bytes();
        let cdr3 = record.get(header_fields["junction"]).unwrap().as_bytes();

        let mut occs = vec![];

        let mut min_dist = 1000;
        'a: for leader in leader_data.iter() {
            let occ: (usize, usize, usize) = find_leading(&leader.myers, &sequence[0..leading_max], args.dist)?;

            if occ.1 > 0 {
                // 亚型一致且错配为0
                if (v_re.is_match(&leader.id)) & (occ.2==0)  {
                    occs.clear();
                    occs.push(occ);
                    _perfect += 1;
                    break 'a;
                }

                // 记录最小错配
                if (min_dist > occ.2) & (occ.1 > 0) {
                    min_dist = occ.2;
                }
                occs.push(occ);
            }
        }

        let (mut fwr1_s, mut fwr1_e):(usize, usize) = (0, 0);
        let (mut fwr2_s, mut fwr2_e):(usize, usize) = (0, 0);
        let (mut fwr3_s, mut fwr3_e):(usize, usize) = (0, 0);
        let (mut fwr4_s, mut fwr4_e):(usize, usize) = (0, 0);

        let mut _leader_start: usize = 0;
        let mut leader_end: usize = 0;
        if !occs.is_empty() {
            found += 1;
            if occs.len() == 1 {
                _leader_start = occs[0].0;
                leader_end = occs[0].1;
            } else {
                let tmp = occs.iter().filter(|&x| x.2==min_dist).collect::<Vec<_>>();
                _leader_start = tmp[0].0;
                leader_end = tmp[0].1;
            }
        }
        if leader_end>0 {
            fwr1_s = leader_end;
        }
        if !cdr1.is_empty() {
            fwr1_e = fwr1_s + sequence.windows(cdr1.len()).skip(fwr1_s).position(|window| window == cdr1).unwrap();
        }
        if fwr1_e>0{
            fwr2_s = fwr1_e + cdr1.len();
        }
        if !cdr2.is_empty(){
            fwr2_e = fwr2_s + sequence.windows(cdr2.len()).skip(fwr2_s).position(|window| window == cdr2).unwrap();
        }

        if fwr2_e>0 {
            fwr3_s = fwr2_e + cdr2.len();
        }
        if !cdr3.is_empty() {
            fwr3_e = fwr3_s + sequence.windows(cdr3.len()).skip(fwr3_s).position(|window| window == cdr3).unwrap();
        }
        if fwr3_e>0 {
            fwr4_s = fwr3_e + cdr3.len();
        }
        if c_data.contains_key(contig_name){
            fwr4_e = *c_data.get(contig_name).unwrap();
        }

        let fwr1 = if (fwr1_s>0) & (fwr1_e>0) {
            from_utf8(&sequence[fwr1_s..fwr1_e]).unwrap()
        } else {
            ""
        };
        let fwr2 = if (fwr2_s>0) & (fwr2_e>0) {
            from_utf8(&sequence[fwr2_s..fwr2_e]).unwrap()
        } else {
            ""
        };
        let fwr3 = if (fwr3_s>0) & (fwr3_e>0) {
            from_utf8(&sequence[fwr3_s..fwr3_e]).unwrap()
        } else {
            ""
        };
        let fwr4 = if (fwr4_s>0) & (fwr4_e>0) {
            from_utf8(&sequence[fwr4_s..fwr4_e]).unwrap()
        } else {
            ""
        };

        wtr.serialize(OutRecord {
            sequence_id: contig_name,
            sequence: from_utf8(sequence).unwrap(),
            sequence_alignment: record.get(header_fields["sequence_alignment"]).unwrap(),
            germline_alignment: record.get(header_fields["germline_alignment"]).unwrap(),
            fwr1,
            cdr1: from_utf8(cdr1).unwrap(),
            fwr2,
            cdr2: from_utf8(cdr2).unwrap(),
            fwr3,
            cdr3: from_utf8(cdr3).unwrap(),
            fwr4,
        })?;
    }
    wtr.flush()?;
    // log::info!("Found {} out of {}; {} perfect.", found, total, perfect);
    log::info!("Found {} out of {}", found, total);
    Ok(())
}
