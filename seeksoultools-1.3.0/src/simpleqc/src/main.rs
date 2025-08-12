use std::thread;
use std::str::FromStr;
use seq_io::fastq::{Reader, OwnedRecord, Record};
use niffler;
use std::collections::{HashSet, HashMap};
use csv::WriterBuilder;
use ndarray::prelude::*;
use ndarray::Array;
use pariter::{IteratorExt as _, scope};
use std::fs::File;
use std::io::{BufWriter, BufReader};
use std::io::prelude::*;
use flate2::Compression;
use flate2::write::GzEncoder;
use std::time::Instant;
use read_structure::ReadStructure;
use anyhow::Result;
use clap::Parser;
use atoi::FromRadix10;

#[macro_use] extern crate log;
extern crate simplelog;

use simplelog::*;


#[derive(Debug)]
pub enum BarcodeState {
    ExactMatch, // 精确匹配
    OneMatch(MutateSeq), // 匹配，错配1
    MultiMatch(Vec<MutateSeq>), // 匹配，错配1
    NoMatch, // 不匹配，
}


/// mutate seq
/// ```
/// let seq = b"A";
/// let bases = b"ACGT";
/// assert_eq!(rfastq::mutate_seq(seq, bases), vec![b"C", b"G", b"T"]);
/// let seq2 = b"AA";
/// let bases = b"ACGT";
/// assert_eq!(rfastq::mutate_seq(seq2, bases), vec![b"CA", b"GA", b"TA", b"AC", b"AG", b"AT"]);
/// ```
pub fn mutate_seq(seq: &[u8], bases: &[u8]) -> Vec<Vec<u8>> {
    let mut v: Vec<Vec<u8>> = Vec::new();
    for i in 0..seq.len() {
        for j in 0..bases.len() {
            if seq[i] != bases[j] {
                let mut new_seq = seq.to_vec();
                new_seq[i] = bases[j];
                v.push(new_seq);
            }
        }
    }
    v
}

/// MutateSeq 结构体
/// 
#[derive(Debug)]
pub struct MutateSeq {
    seq: Vec<u8>,
    pos: usize,
    raw_base: u8,
    new_base: u8
}

pub fn mutate_seq2(seq: &[u8], bases: &[u8]) -> Vec<MutateSeq> {
    let mut v: Vec<MutateSeq> = Vec::new();
    for i in 0..seq.len() {
        for j in 0..bases.len() {
            if seq[i] != bases[j] {
                let mut new_seq = seq.to_vec();
                new_seq[i] = bases[j];
                v.push(MutateSeq{
                    seq: new_seq,
                    pos: i,
                    raw_base: seq[i],
                    new_base: bases[j]
                });
            }
        }
    }
    v
}
pub fn search_whitelist(barcode: &[u8], white_list: &HashSet<Vec<u8>>, mismatch_allowed:bool) -> BarcodeState{
    if white_list.contains(barcode) {
        return BarcodeState::ExactMatch;
    }
    if mismatch_allowed {
        let barcodes = mutate_seq2(barcode, b"ACGT");
        let mut barcode_mismatch: Vec<MutateSeq> = vec![];
        for s in barcodes {
            if white_list.contains(&s.seq) {
                // return BarcodeState::Match(i);
                barcode_mismatch.push(s)
            }
        }
        let candidates_num = barcode_mismatch.len();
        return  match candidates_num  {
            0 => BarcodeState::NoMatch,
            1 => BarcodeState::OneMatch(barcode_mismatch.into_iter().nth(0).unwrap()),
            _ => BarcodeState::MultiMatch(barcode_mismatch),
        }
    }
    return BarcodeState::NoMatch;
}
fn read_white_list(f: &str)-> Result<HashSet<Vec<u8>>> {
    let mut white_list: HashSet<Vec<u8>> = HashSet::new();
    // let file = File::open(f).unwrap();
    let (reader, _format) = niffler::send::from_path(f)?;
    let reader = BufReader::new(reader);

    for line in reader.lines() {
        let barcode = line.unwrap().as_bytes().to_owned();
        white_list.insert(barcode);
    }
    Ok(white_list)
}

// 定义一个用于存储Read状态的结构体
#[derive(Debug)]
struct ReadStat{
    // 存储序列相关的统计信息，采用二维数组形式，其中的具体类型和维度大小在初始化时确定
    seq_stat: Array<usize, Dim<[usize; 2]>>,
    // 存储质量相关的统计信息，同样采用二维数组形式，允许与序列统计信息不同的维度大小
    qual_stat: Array<usize, Dim<[usize; 2]>>,
}

// 实现ReadStat结构体的操作
impl ReadStat {
    /// 创建一个新的ReadStat实例。
    ///
    /// 初始化序列统计（seq_stat）和质量统计（qual_stat）数组，为后续的统计操作准备。
    fn new() -> Self {
        let seq_stat = Array::<usize, _>::zeros((512, 5).f());
        let qual_stat = Array::<usize, _>::zeros((512, 50).f());
        ReadStat {
            seq_stat,
            qual_stat,
        }
    }

    /// 更新一个ReadStat实例的统计信息。
    ///
    /// 将两个ReadStat实例的序列统计（seq_stat）和质量统计（qual_stat）相加，
    /// 结果存储在第一个实例（a）中。
    fn update(a: &mut ReadStat, b:ReadStat){
        a.seq_stat = &a.seq_stat + b.seq_stat;
        a.qual_stat = &a.qual_stat + b.qual_stat;
    }
}
/// `Chunked` 结构体用于对迭代器进行分块处理。
/// 
/// 分块处理是一种将大量数据分割成更小、更易于管理的块的操作，
/// 这对于处理大数据集时非常有用，例如读取大文件或网络传输中的数据。
/// 
/// # 参数
/// - `I`: 一个可以转换为迭代器的类型。这允许`Chunked`处理各种可以迭代的数据结构。
/// - `chunk_size`: 分块的大小，决定了每块包含的元素数量。如果一个块的元素数量不足`chunk_size`，
///   那么这是最后一个块，之后不会再有块。
/// 
/// # 示例
/// 
/// ```rust
/// use std::iter::IntoIterator;
/// 
/// let data = vec![1, 2, 3, 4, 5];
/// let chunked = Chunked {
///     iterator: data.into_iter(),
///     chunk_size: 2,
/// };
/// 
/// for chunk in chunked.collect::<Vec<_>>() {
///     println!("{:?}", chunk);
/// }
/// ```
/// 
/// 这将输出:
/// ```
/// [1, 2]
/// [3, 4]
/// [5]
/// ```
struct Chunked<I>
where I: IntoIterator
{
    iterator: I,
    chunk_size: usize,
}
/// 创建一个分块迭代器，将集合按照指定大小进行分块
///
/// # 参数
/// - `c`: 一个可转换为迭代器的集合
/// - `chunk_size`: 每个分块的大小
///
/// # 返回
/// 返回一个`Chunked`结构体，其中包含原始迭代器和分块大小
///
/// # 类型参数
/// `Collection`: 一个实现`IntoIterator` trait 的集合类型
///
/// # 类型约束
/// - `Collection: IntoIterator`: 集合类型必须能够转换为迭代器
///
/// # 例子
/// ```rust
/// let numbers = vec![1, 2, 3, 4, 5, 6, 7, 8];
/// let chunks = chunked(numbers, 3);
/// // chunks现在是一个分块迭代器，每个分块大小为3
/// ```
fn chunked<Collection>(c: Collection, chunk_size: usize) -> Chunked<Collection::IntoIter>
where
    Collection: IntoIterator,
{
    // 获取集合的迭代器
    let iterator = c.into_iter();
    
    // 创建并返回一个Chunked结构体，设置迭代器和分块大小
    Chunked {
        iterator,
        chunk_size,
    }
}

/// 实现迭代器的分块功能
///
/// # 泛型参数
/// - `I`: 一个实现了迭代器 trait 的类型
///
/// # 类型关联
/// - `Item`: 一个包装了原始迭代器元素类型向量的类型
///
/// # 方法
/// - `next`: 返回一个 Option 类型，其中包含最多包含 `chunk_size` 个元素的向量
///   - 该方法通过从底层迭代器中按需获取元素来创建分块
///   - 每次调用都会从迭代器中获取最多 `chunk_size` 个元素，并将它们收集到一个向量中
///   - 如果获取的元素为空向量，则过滤掉，不返回
impl<I: Iterator> Iterator for Chunked<I> {
    type Item = Vec<I::Item>;

    fn next(&mut self) -> Option<Self::Item> {
        Some(self.iterator.by_ref().take(self.chunk_size).collect())
            .filter(|chunk: &Vec<_>| !chunk.is_empty())
    }
}
/// 将碱基转换为对应的数字
///
/// # 参数
/// - `base`: 一个代表碱基的 u8 类型的引用
///
/// # 返回值
/// - `Some(u8)`：表示成功转换的数字，范围为 0 到 4，分别对应不同的碱基
/// - `None`：如果输入的碱基无法识别，则返回 None
///
/// # 说明
/// 本函数用于将 DNA/RNA 中的碱基符号转换为数字编码，以方便计算机处理和存储
/// 特别地，'T' 和 'U' 被视为相同的碱基，因为在 RNA 中，'T' 会被 'U' 取代
/// 'N' 代表任何未知或未指定的碱基，而 None 则用于表示无法识别的输入
fn base_to_number(base: &u8) -> Option<u8> {
    match base {
        &b'A' => Some(0),
        &b'T' | &b'U' => Some(1), // U 代表 Uracil，在 RNA 中使用
        &b'C' => Some(2),
        &b'G' => Some(3),
        &b'N' => Some(4), // N 代表未知或任意碱基
        _ => None,       // 如果不是上述碱基，则返回 None
    }
}
fn write_read(read:OwnedRecord, header: &Vec<u8>, buffer:&mut Vec<u8>){
    buffer.write_all(b"@").unwrap();
    let mut new_header = header.clone();
    new_header.extend_from_slice(b"_");
    new_header.extend_from_slice(&read.head());
    buffer.write_all(&new_header).unwrap();
    buffer.write_all(b"\n").unwrap();
    buffer.write_all(read.seq()).unwrap();
    buffer.write_all(b"\n+\n").unwrap();
    buffer.write_all(read.qual()).unwrap();
    buffer.write_all(b"\n").unwrap();
}

fn compress_data(buffer: &Vec<u8>)-> Vec<u8>{
    let mut e = GzEncoder::new(Vec::new(), Compression::default());
    e.write_all(buffer).unwrap();
    e.finish().unwrap()
}

fn read_stat(
    chunk: Vec<(Result<OwnedRecord, seq_io::fastq::Error>, Result<OwnedRecord, seq_io::fastq::Error>)>,
    rs: &ReadStructure,
    wl: &HashSet<Vec<u8>>
// ) -> ((ReadStat, ReadStat), (usize, usize))
) -> ((ReadStat, ReadStat), (Vec<u8>, Vec<u8>), (usize, usize), (Vec<u8>, Vec<u8>), HashMap<Vec<u8>, usize>)
{
    let mut r1_stat = ReadStat::new();
    let mut r2_stat = ReadStat::new();
    let mut barcode_stat: HashMap<Vec<u8>, usize> = HashMap::new();

    let mut total_reads_num: usize= 0;
    let mut valid_reads_num: usize = 0;
    let mut buffer1: Vec<_> =vec![];
    let mut buffer2: Vec<_> = vec![];
    let mut multi_buffer1: Vec<_> =vec![];
    let mut multi_buffer2: Vec<_> = vec![];
    for (r1, r2) in chunk {
        let r1 = r1.unwrap();
        let r2 = r2.unwrap();
        for i in 0..r1.seq.len() {
            let base = r1.seq[i];
            let base_num = base_to_number(&base).unwrap() as usize;
            let qual = r1.qual[i] as usize - 33;
            r1_stat.seq_stat[[i, base_num]] += 1;
            r1_stat.qual_stat[[i, qual]] += 1;
        }
        for i in 0..r2.seq.len() {
            let base = r2.seq[i];
            let base_num = base_to_number(&base).unwrap() as usize;
            let qual = r2.qual[i] as usize - 33;
            r2_stat.seq_stat[[i, base_num]] += 1;
            r2_stat.qual_stat[[i, qual]] += 1;
        }
        total_reads_num += 1;
        let cb = rs.cellular_barcodes().next().unwrap()
            .extract_bases(&r1.seq()).unwrap();
        match search_whitelist(cb, wl, true){
            BarcodeState::ExactMatch => {
                let mut  new_header = Vec::new();
                let umi = rs.molecular_barcodes().next().unwrap()
                    .extract_bases(&r1.seq()).unwrap();
                new_header.extend_from_slice(cb);
                new_header.extend_from_slice(b"_");
                new_header.extend_from_slice(umi);
                new_header.extend_from_slice(b"_M");
                valid_reads_num += 1;
                barcode_stat.entry(cb.to_vec()).and_modify(|n| *n += 1).or_insert(1);
                write_read(r1, &new_header,&mut buffer1);
                write_read(r2, &new_header,&mut buffer2);
                
            },
            BarcodeState::OneMatch(m) => {
                let mut new_header = Vec::new();
                let umi = rs.molecular_barcodes().next().unwrap()
                    .extract_bases(&r1.seq()).unwrap();
                valid_reads_num += 1;
                new_header.extend_from_slice(&m.seq);
                new_header.extend_from_slice(b"_");
                new_header.extend_from_slice(umi);
                new_header.extend_from_slice(b"_");
                new_header.extend_from_slice(m.pos.to_string().as_bytes());
                new_header.push(m.raw_base);
                barcode_stat.entry(m.seq).and_modify(|n| *n += 1).or_insert(1);
                write_read(r1, &new_header,&mut buffer1);
                write_read(r2, &new_header,&mut buffer2);    
            },
            BarcodeState::MultiMatch(m) => {
                let mut new_header = Vec::new();
                let umi = rs.molecular_barcodes().next().unwrap()
                    .extract_bases(&r1.seq()).unwrap();
                // valid_reads_num += 1;
                new_header.extend_from_slice(cb);
                new_header.extend_from_slice(b"_");
                let mut a = m.into_iter().fold(Vec::new(), |mut acc, m| {
                    let pos = m.pos.to_string().as_bytes().to_vec();
                    acc.extend_from_slice(&pos);
                    acc.push(m.new_base);
                    acc.extend_from_slice(b":");
                    acc
                });
                a.pop();
                new_header.extend_from_slice(&a);
                new_header.extend_from_slice(b"_");
                new_header.extend_from_slice(umi);
                // let mut new_header1 = new_header.clone();
                // let mut new_header2 = new_header;
                // new_header1.extend_from_slice(r1.head());
                // new_header2.extend_from_slice(r2.head());
                write_read(r1, &new_header,&mut multi_buffer1);
                write_read(r2, &new_header,&mut multi_buffer2);
            },
            BarcodeState::NoMatch => {},
        }
    }
    // ((r1_stat, r2_stat), (valid_reads_num, total_reads_num))
    (
        (r1_stat, r2_stat),
        (compress_data(&buffer1), compress_data(&buffer2)),
        (valid_reads_num, total_reads_num),
        (compress_data(&multi_buffer1), compress_data(&multi_buffer2)),
        barcode_stat
    )
}

fn process_multi(
    chunk: Vec<(Result<OwnedRecord, seq_io::fastq::Error>, Result<OwnedRecord, seq_io::fastq::Error>)>,
    barcode_stat: &HashMap<Vec<u8>, usize>
)->((Vec<u8>, Vec<u8>), (usize, usize)){

    let mut total_reads_num: usize = 0;
    let mut valid_reads_num: usize = 0;
    let mut buffer1: Vec<_> =vec![];
    let mut buffer2: Vec<_> = vec![];
    for (r1, r2) in chunk {
        let r1 = r1.unwrap();
        let r2 = r2.unwrap();
        total_reads_num += 1;
        let head1 = r1.head().split(|b| *b==b'_').collect::<Vec<&[u8]>>();
        let (raw_barcode, info, umi) = (head1[0], head1[1], head1[2]);

        let mut new_barcodes = info.split(|b| *b==b':').collect::<Vec<&[u8]>>()
            .into_iter().map(|x| {
                let (num, idx) =  usize::from_radix_10(x);
                let mut new_barcode = raw_barcode.to_vec();
                new_barcode[num] = x[idx];
                let count = match barcode_stat.get(&new_barcode){
                    Some(x) => *x,
                    None => 0
                };
                (new_barcode, count, x)
            }).collect::<Vec<(Vec<u8>, usize, &[u8])>>();
        new_barcodes
            .sort_by(|a, b| b.1.cmp(&a.1));

        if new_barcodes[0].1 > new_barcodes[1].1 {
            valid_reads_num += 1;
            let mut new_header =  Vec::new();
            new_header.extend(new_barcodes[0].0.clone());
            new_header.extend_from_slice(b"_");
            new_header.extend_from_slice(umi);
            new_header.extend_from_slice(b"_");
            new_header.extend_from_slice(new_barcodes[0].2);

            // new_header1 = new_header.clone();
            // new_header2 = new_header;
            // new_header1.extend_from_slice(last1);
            // new_header2.extend_from_slice(last2);
            write_read(r1, &new_header,&mut buffer1);
            write_read(r2, &new_header,&mut buffer2);   
        }
    }
    (
        (compress_data(&buffer1), compress_data(&buffer2)),
        (valid_reads_num, total_reads_num),
    )
}

/// 将序列统计信息和质量统计信息写入到CSV文件中
/// 
/// # 参数
/// - `stat`: 借用的读取统计信息，包含序列和质量统计数据
/// - `out_prefix`: 输出文件名前缀，用于生成输出文件名
fn write_stat(stat: &ReadStat, out_prefix: String) {
    // 构造序列统计信息的文件名，并创建相应的CSV写入器
    let seq_stat_file = out_prefix.clone() + "base_stat.csv";
    let mut wtr0 = WriterBuilder::new().from_path(seq_stat_file).unwrap();
    // 写入序列统计信息的表头
    wtr0.write_record(&["pos", "A", "T", "C", "G", "N"]).unwrap();
    
    // 找到序列统计信息中有非零值的第一行索引
    let mut row_idx: usize = 0;
    for (i, r) in stat.seq_stat.rows().into_iter().enumerate() {
        if r.into_iter().any(|x| x>&0) {
            row_idx = i;
        }
    }
    
    // 定义列索引，并切片写入序列统计信息
    let col_idx: usize = 4;
    stat.seq_stat.slice(s![0..row_idx+1, 0..col_idx+1]).rows().into_iter().enumerate().for_each(|(idx, row)| {
        let mut row_data = vec![idx.to_string()];
        let mut data: Vec<String> = row.iter().map(|&v| v.to_string()).collect();
        row_data.append(&mut data);
        wtr0.write_record(&row_data).unwrap();
    });
    
    // 构造质量统计信息的文件名，并创建相应的CSV写入器
    let qual_stat_file = out_prefix.clone() + "qual_stat.csv";
    let mut wtr1 = WriterBuilder::new().from_path(qual_stat_file).unwrap();
    
    // 找到质量统计信息中有非零值的第一行和第一列索引
    let mut row_idx: usize = 0;
    for (i, r) in stat.qual_stat.rows().into_iter().enumerate() {
        if r.into_iter().any(|x| x>&0) {
            row_idx = i;
        }
    }
    
    let mut col_idx: usize = 0;
    for (j, c) in stat.qual_stat.columns().into_iter().enumerate() {
        if c.into_iter().any(|x| x>&0) {
            col_idx = j;
        }
    }
    
    // 写入质量统计信息的表头
    let mut header = vec!["pos".to_string()];
    let mut qual = (0..col_idx + 1).map(|x| x.to_string()).collect::<Vec<String>>();
    header.append(&mut qual);
    wtr1.write_record(&header).unwrap();
    
    // 切片并写入质量统计信息
    stat.qual_stat.slice(s![0..row_idx+1, 0..col_idx+1]).rows().into_iter().enumerate().for_each(|(idx, row)| {
        let mut row_data = vec![idx.to_string()];
        let mut data: Vec<String> = row.iter().map(|&v| v.to_string()).collect();
        row_data.append(&mut data);
        wtr1.write_record(&row_data).unwrap();
    });
    
    // 刷新写入器，确保所有数据都写入到文件中
    wtr1.flush().unwrap();
}


/// Simple program to stat fq
#[derive(Parser, Debug)]
#[command(version, about, long_about = None, arg_required_else_help(true))]
struct Args {
    /// read1 fastq file
    #[arg(long, required(true), num_args = 1..)]
    fq1: Vec<String>,

    /// read2 fastq file
    #[arg(long, required(true), num_args = 1..)]
    fq2: Vec<String>,

    /// read1 structure
    #[arg(long, required(true))]
    rs: String,

    /// white list file
    #[arg(short, long, required(true))]
    wl: String,
    
    /// Chunk size
    #[arg(short, long, default_value_t = 500000)]
    chunk: usize,

    /// Threads num
    #[arg(short, long, default_value_t = 4)]
    threads: usize,

    /// Output file prefix
    #[arg(short, long, default_value_t = String::from("./stat_"))]
    out: String,
}

fn main() -> Result<()> {
    let _ = SimpleLogger::init(
        LevelFilter::Info,
        ConfigBuilder::new()
            .set_time_format_rfc3339()
            .set_time_offset_to_local().unwrap()
            .build()
    );
    let now = Instant::now();
    let args = Args::parse();
    let count = thread::available_parallelism()?.get();
    let mut threads: usize = args.threads;
    if threads > count {
        info!("Threads num is too large, set to {}", count);
        threads = count;
    }
    let workers = if threads >1 {threads} else {1};
    info!("available parallelism: {count}, workers: {workers}");
    let paths = &args.fq1.into_iter().zip(&args.fq2).collect::<Vec<_>>();

    let white_list = read_white_list(&args.wl)?;

    let rs = ReadStructure::from_str(&args.rs)?;

    let r1_data = & mut ReadStat::new();
    let r2_data = & mut ReadStat::new();

    let mut total_num: usize = 0;
    let mut valid_num: usize = 0;
    let mut barcode_stat: HashMap<Vec<u8>, usize> = HashMap::new();

    let f1 = File::create(args.out.clone() + "r1.fq.gz")?;
    let mut w1 = BufWriter::new(f1);
    let f2 = File::create(args.out.clone() + "r2.fq.gz")?;
    let mut w2 = BufWriter::new(f2);

    let multi_f1 = File::create(args.out.clone() + "r1.multi.fq.gz")?;
    let mut multi_w1 = BufWriter::new(multi_f1);
    let multi_f2 = File::create(args.out.clone() + "r2.multi.fq.gz")?;
    let mut multi_w2 = BufWriter::new(multi_f2);
    for (fq1, fq2) in paths {
        info!("Processing file: {}, {}", fq1, fq2);
        let (reader1, _format) = niffler::send::from_path(fq1)?;
        let (reader2, _format) = niffler::send::from_path(fq2)?;

        let reader1 = Reader::new(reader1);
        let reader2 = Reader::new(reader2);

        let reader = reader1.into_records().zip(reader2.into_records());
        let _ = scope(|scope| {
            chunked(reader, args.chunk)
                .parallel_map_scoped_custom(scope,|o| o.threads(workers), |chunk|{
                    read_stat(chunk, &rs, &white_list)
                })
                .for_each(|(
                    (r1_stat, r2_stat),
                    (buf1, buf2),
                    (valid_reads_num, total_reads_num),
                    (multi_buf1, multi_buf2),
                    chunk_barcode_stat

                )| {
                // .for_each(|((r1_stat, r2_stat), (valid_reads_num, total_reads_num))| {
                    ReadStat::update(r1_data, r1_stat);
                    ReadStat::update(r2_data, r2_stat);
                    valid_num += valid_reads_num;
                    total_num += total_reads_num;
                    w1.write_all(&buf1).unwrap();
                    w2.write_all(&buf2).unwrap();
                    multi_w1.write_all(&multi_buf1).unwrap();
                    multi_w2.write_all(&multi_buf2).unwrap();
                    chunk_barcode_stat.iter().for_each(|(k, v)| {
                        barcode_stat.entry(k.to_vec()).and_modify(|n| *n += v).or_insert(*v);
                    });
                    // info!("{} reads processed, {} reads valid.", total_num, valid_num);
                })
        }); 
    }
    multi_w1.flush()?;
    multi_w2.flush()?;
    info!("Processing multiple candidates...");

    let (reader1, _format) = niffler::send::from_path(args.out.clone() + "r1.multi.fq.gz")?;
    let (reader2, _format) = niffler::send::from_path(args.out.clone() + "r2.multi.fq.gz")?;
    let reader1 = Reader::new(reader1);
    let reader2 = Reader::new(reader2);
    let reader = reader1.into_records().zip(reader2.into_records());
    let _ = scope(|scope| {
        chunked(reader, args.chunk)
            .parallel_map_scoped_custom(scope,|o| o.threads(workers), |chunk|{
                process_multi(chunk, &barcode_stat)
            })
            .for_each(|(
                (buf1, buf2),
                (valid_reads_num, _total_reads_num),
            )| {
                w1.write_all(&buf1).unwrap();
                w2.write_all(&buf2).unwrap();
                valid_num += valid_reads_num;
            })
    });

    w1.flush()?;
    w2.flush()?;

    // println!("{:?}", barcode_stat);
    write_stat(&r1_data, args.out.clone() + "r1_");
    write_stat(&r2_data, args.out.clone() + "r2_");
    let data = format!("total_num: {total_num}\nvalid_num: {valid_num}\n");
    File::create(args.out.clone() + "summary.txt")?.write_all(data.as_bytes())?;
    info!("Elapsed time: {} seconds", now.elapsed().as_secs());
    Ok(())
}
