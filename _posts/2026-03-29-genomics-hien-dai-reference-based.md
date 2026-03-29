---
layout: post
title:  "Tư Duy Lớn Của Genomics Hiện Đại: Reference Genome Và Hệ Sinh Thái Công Nghệ Đa Omics"
author: dat
categories: [ Bioinformatics, Genomics ]
image: assets/my_figs/ds/genomics_hien_dai.png
tags: [ genomics, rnaseq, epigenomics, multi-omics, sequencing, reference-genome, bioinformatics ]
---

Nếu bạn mới bắt đầu với tin sinh học, bạn sẽ sớm nhận ra rằng gần như mọi phân tích đều xoay quanh một thứ duy nhất: **bộ gen tham chiếu** (reference genome). Đây không chỉ là một file dữ liệu — đây là hệ tọa độ của toàn bộ sinh học phân tử hiện đại. Hiểu được tư duy này sẽ giúp bạn thấy được bức tranh toàn cảnh: tại sao cần nhiều công nghệ khác nhau, chúng kết nối với nhau như thế nào, và làm sao chúng ta có thể nghiên cứu hàng chục "lớp" thông tin sinh học chỉ từ một mẫu tế bào.

## 1. Reference Genome — Hệ Tọa Độ Của Sinh Học Phân Tử

Hãy tưởng tượng bộ gen là một cuốn sách khổng lồ — ~3 tỷ chữ cái (base pair) đối với người. Mỗi chữ cái là một nucleotide (A, T, G, C). Dự án Human Genome Project hoàn thành năm 2003 lần đầu tiên "đọc" được toàn bộ cuốn sách này, tạo ra **reference genome** — bản tham chiếu chuẩn cho loài người (GRCh38/hg38).

**Tại sao reference genome lại quan trọng đến vậy?**

Vì nó cho phép chúng ta **so sánh**. Thay vì nghiên cứu mỗi cá thể từ đầu, ta chỉ cần hỏi: *"Bộ gen của người bệnh này khác với bản tham chiếu ở chỗ nào? Gen nào đang được bật/tắt hơn mức bình thường? Vùng DNA nào đang mở hay đóng trong tế bào ung thư?"*

Mọi công nghệ omics hiện đại — dù là giải trình tự DNA, RNA, hay đo methylation — đều hoạt động theo một quy trình chung:

```
Thu thập mẫu sinh học
       ↓
Phân tách phân tử quan tâm (DNA / RNA / protein...)
       ↓
Chuyển đổi thành tín hiệu đọc được (thường là trình tự DNA)
       ↓
"Map" (căn chỉnh) kết quả lên reference genome
       ↓
Đếm, so sánh, phân tích thống kê
```

Bước **"map lên reference"** chính là bước biến tất cả các công nghệ rời rạc thành một ngôn ngữ chung — ngôn ngữ của tọa độ nhiễm sắc thể.

---

## 2. Công Nghệ Giải Trình Tự — Nền Tảng Của Tất Cả

Để đọc được trình tự DNA/RNA, ta cần **công nghệ giải trình tự (sequencing)**. Có hai thế hệ chính đang được dùng song song:

### Short-read Sequencing (Illumina)

Illumina là công nghệ phổ biến nhất hiện nay. Nguyên lý là **sequencing by synthesis**: DNA được nhân bản thành hàng triệu cụm (cluster) trên một chip, rồi đọc từng nucleotide thông qua tín hiệu huỳnh quang.

- **Độ dài đọc**: 100–300 bp mỗi read
- **Thông lượng**: rất cao, hàng chục Gb mỗi lane
- **Độ chính xác**: ~99.9% (Q30)
- **Dùng cho**: WGS, RNA-seq, ATAC-seq, ChIP-seq, WGBS

```
Ưu điểm: Rẻ, nhanh, chính xác cao
Nhược điểm: Read ngắn → khó assembly vùng lặp lại, khó phát hiện biến thể cấu trúc lớn
```

### Long-read Sequencing (PacBio & Oxford Nanopore)

Thế hệ thứ ba đọc từng phân tử DNA đơn mà không cần nhân bản:

- **PacBio HiFi (SMRT)**: read dài 10–25 kb, độ chính xác ~99.9%
- **Oxford Nanopore (ONT)**: read dài đến hàng Mb (toàn nhiễm sắc thể!), real-time, portable

```
Ưu điểm: Đọc qua vùng lặp lại, phát hiện methylation trực tiếp, assembly T2T
Nhược điểm: Đắt hơn (PacBio), hoặc tỉ lệ lỗi cao hơn và cần điều chỉnh (Nanopore raw)
```

**Kết hợp cả hai** (hybrid assembly) là chiến lược vàng cho genome assembly: long-read cung cấp khung xương, short-read "đánh bóng" (polish) độ chính xác.

---

## 3. Các Lớp Omics — Đọc Nhiều Chiều Từ Cùng Một Hệ Tọa Độ

Đây là điểm mấu chốt: từ một tế bào, có vô số lớp thông tin sinh học. Mỗi công nghệ "soi" vào một lớp khác nhau, nhưng tất cả đều được quy về cùng một tọa độ trên reference genome.

### 🧬 Genomics — Lớp DNA

**Mục tiêu**: Tìm biến thể di truyền (SNP, indel, CNV, SV) trong DNA hệ gen.

**Công nghệ**: Whole Genome Sequencing (WGS), Whole Exome Sequencing (WES)

**Pipeline chuẩn** (GATK Best Practices):
```bash
# 1. Căn chỉnh lên reference
bwa mem -t 8 hg38.fa sample_R1.fastq.gz sample_R2.fastq.gz | \
    samtools sort -o sample.sorted.bam

# 2. Đánh dấu duplicate reads
gatk MarkDuplicates -I sample.sorted.bam -O sample.dedup.bam \
    -M sample.dedup_metrics.txt

# 3. Base Quality Score Recalibration
gatk BaseRecalibrator -I sample.dedup.bam -R hg38.fa \
    --known-sites dbsnp_146.hg38.vcf.gz -O recal.table

# 4. Gọi biến thể
gatk HaplotypeCaller -I sample.dedup.bam -R hg38.fa \
    -O sample.g.vcf.gz -ERC GVCF
```

**Câu hỏi trả lời được**: Người bệnh này có đột biến BRCA1 không? Gen này có bị mất đoạn (deletion) trong tế bào ung thư không?

---

### 🔬 Transcriptomics — Lớp RNA

**Mục tiêu**: Đo lường những gen nào đang được biểu hiện (express) và ở mức độ nào.

**Công nghệ**: RNA-seq (bulk), Single-cell RNA-seq (scRNA-seq)

RNA không thể giải trình tự trực tiếp — nó được chuyển ngược thành cDNA (complementary DNA) trước. Sau đó cDNA được giải trình tự bằng Illumina và map lên reference transcriptome.

**Pipeline RNA-seq đơn giản**:
```bash
# 1. Kiểm tra chất lượng
fastqc sample_R1.fastq.gz sample_R2.fastq.gz

# 2. Căn chỉnh siêu nhanh với STAR
STAR --runThreadN 8 \
     --genomeDir /ref/hg38_star \
     --readFilesIn sample_R1.fastq.gz sample_R2.fastq.gz \
     --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix sample_

# 3. Đếm reads theo gen
featureCounts -T 8 -a hg38.gtf -o counts.txt sample_Aligned.sortedByCoord.out.bam
```

**Phân tích biểu hiện gen vi sai (DEG) với R/DESeq2**:
```r
library(DESeq2)

# Load count matrix
counts <- read.table("counts.txt", header=TRUE, row.names=1)
coldata <- data.frame(condition = c("control","control","tumor","tumor"))

# Tạo DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts,
                               colData   = coldata,
                               design    = ~ condition)

# Chạy phân tích
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "tumor", "control"))

# Lọc gen significant
sig_genes <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
```

**scRNA-seq** — đưa mọi thứ lên tầm cao mới: thay vì trung bình hoá toàn bộ mô, ta đọc biểu hiện gen của **từng tế bào đơn lẻ** (hàng nghìn đến hàng triệu tế bào). Điều này giúp phát hiện các kiểu tế bào (cell types), trạng thái tế bào, và quỹ đạo phát triển (trajectory).

---

### 🔓 Epigenomics — Lớp "Bật/Tắt" Gen

DNA không tồn tại trần trụi — nó cuộn quanh protein histone, tạo thành cấu trúc **chromatin**. Trạng thái "mở" hay "đóng" của chromatin quyết định gen nào được phép biểu hiện. Đây là lớp **epigenomics**.

**Công nghệ chính**:

| Công nghệ | Đo lường |
|-----------|----------|
| **ATAC-seq** | Vùng chromatin mở (accessibility) |
| **ChIP-seq** | Liên kết protein/histone với DNA |
| **WGBS / RRBS** | Methylation DNA (5-methylcytosine) |
| **CUT&RUN / CUT&TAG** | Histone marks, transcription factor binding |
| **Hi-C** | Cấu trúc 3D nhiễm sắc thể |

**ATAC-seq** là một ví dụ tiêu biểu: enzyme Tn5 transposase chỉ có thể "cắt" vào những vùng DNA đang mở (không bị nucleosome che). Các đoạn DNA bị cắt được thu thập và giải trình tự — vùng nào có nhiều reads chính là vùng chromatin đang mở, tức vùng **điều hòa gene đang hoạt động**.

```bash
# ATAC-seq pipeline cơ bản
# 1. Trim adapters
trimmomatic PE sample_R1.fastq.gz sample_R2.fastq.gz \
    trimmed_R1.fastq.gz /dev/null trimmed_R2.fastq.gz /dev/null \
    ILLUMINACLIP:adapters.fa:2:30:10

# 2. Căn chỉnh (loại bỏ mitochondrial reads)
bowtie2 -X 2000 -x hg38 -1 trimmed_R1.fastq.gz -2 trimmed_R2.fastq.gz | \
    samtools view -bS -F 4 | samtools sort -o atac.bam

# 3. Gọi đỉnh (peak calling)
macs2 callpeak -t atac.bam -f BAMPE -n sample \
    --nomodel --shift -100 --extsize 200 \
    --outdir peaks/
```

---

### 🧫 Multi-omics — Khi Tất Cả Hội Tụ

Sức mạnh thực sự của genomics hiện đại nằm ở **tích hợp đa lớp omics** từ cùng một mẫu (hoặc đối tượng nghiên cứu). Mỗi lớp trả lời một câu hỏi khác nhau:

```
┌─────────────────────────────────────────────────────┐
│                 Reference Genome (hg38)             │
│              (Hệ tọa độ dùng chung)                 │
└──────────┬────────┬────────┬────────┬───────────────┘
           │        │        │        │
     Genomics  Transcr.  Epigen.  Proteomics
     (WGS)    (RNA-seq) (ATAC-seq)(Mass Spec)
           │        │        │        │
       Các biến  Biểu hiện  Vùng    Protein
       thể di    gen        mở      thực sự
       truyền               chromatin có mặt
           │        │        │        │
           └────────┴────────┴────────┘
                          │
                  Câu hỏi sinh học:
              "Tại sao tế bào này trở thành
               tế bào ung thư?"
```

**Ví dụ tích hợp thực tế**: Trong nghiên cứu ung thư vú:
1. **WGS** → phát hiện đột biến somatic trong PIK3CA
2. **RNA-seq** → xác nhận con đường PI3K/AKT bị upregulate
3. **ATAC-seq** → thấy vùng promoter của gen kháng thuốc đang mở
4. **ChIP-seq H3K27ac** → xác nhận enhancer đang active
5. **scRNA-seq** → phát hiện subpopulation tế bào kháng thuốc chiếm 5% khối u

Kết hợp 5 lớp thông tin này, nhà nghiên cứu có thể đề xuất cơ chế kháng thuốc và thiết kế chiến lược điều trị chính xác.

---

## 4. Công Cụ Bioinformatics — Hệ Sinh Thái Pipeline

Mỗi bước phân tích có những công cụ tiêu chuẩn được cộng đồng tin tưởng:

| Bước | Công cụ phổ biến |
|------|------------------|
| Kiểm tra QC | FastQC, MultiQC |
| Trim adapter | Trimmomatic, fastp, Cutadapt |
| Short-read alignment | BWA-MEM2, STAR (RNA), Bowtie2 |
| Long-read alignment | Minimap2 |
| BAM processing | SAMtools, Picard |
| Variant calling | GATK, DeepVariant, FreeBayes |
| DE analysis | DESeq2, edgeR (R) |
| Peak calling | MACS2, MACS3 |
| scRNA-seq | Seurat (R), Scanpy (Python) |
| Multi-omics | Seurat WNN, MOFA+, mixOmics |
| Workflow manager | Snakemake, Nextflow, nf-core |

**Nextflow/nf-core** đặc biệt đáng chú ý: đây là các pipeline chuẩn hoá, đã được kiểm thử, dễ chạy trên cả máy tính cá nhân lẫn HPC hay cloud:

```bash
# Chạy pipeline RNA-seq chuẩn nf-core với một lệnh
nextflow run nf-core/rnaseq \
    --input samplesheet.csv \
    --outdir results/ \
    --genome GRCh38 \
    -profile docker
```

---

## 5. Tư Duy Tổng Thể — Từ Câu Hỏi Sinh Học Đến Insight

Điều quan trọng cần nhớ không phải là thuộc lòng từng công nghệ, mà là **tư duy ngược từ câu hỏi sinh học**:

1. **Câu hỏi là gì?** → Tế bào ung thư này dùng con đường nào để sống sót?
2. **Cần lớp thông tin nào?** → Biến thể DNA? Biểu hiện gen? Epigenetics?
3. **Công nghệ nào phù hợp?** → WGS + RNA-seq + ATAC-seq?
4. **Pipeline nào chuẩn nhất?** → nf-core/sarek (somatic), nf-core/rnaseq, nf-core/atacseq
5. **Tích hợp dữ liệu như thế nào?** → MOFA+, Seurat WNN, custom R/Python scripts

Reference genome không phải là điểm bắt đầu của một phân tích — nó là **cái nền** mà tất cả mọi thứ được dán lên. Mỗi bead trong chuỗi hạt multi-omics đều biết vị trí của mình vì tất cả đều dùng chung một hệ tọa độ: **nhiễm sắc thể:vị trí:sợi** (chr:pos:strand).

---

## Kết Luận

Genomics hiện đại không phải là một công nghệ duy nhất — đó là một **hệ sinh thái** gồm nhiều công nghệ đọc các lớp thông tin khác nhau, tất cả được neo vào reference genome như một tọa độ chung. Hiểu được mối liên hệ này giúp bạn:

- Chọn đúng công nghệ cho câu hỏi sinh học của mình
- Hiểu tại sao các pipeline bioinformatics được thiết kế theo cách đó
- Tích hợp nhiều loại dữ liệu để trả lời câu hỏi phức tạp hơn

Trong các bài viết tiếp theo, chúng ta sẽ đi sâu vào từng lớp omics: bắt đầu với RNA-seq từ A đến Z, rồi đến ATAC-seq, và cuối cùng là single-cell multi-omics.

**Tài liệu tham khảo:**
- [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows){:target="_blank"}
- [nf-core pipelines](https://nf-co.re/pipelines){:target="_blank"}
- [Bioconductor — R packages cho bioinformatics](https://bioconductor.org/){:target="_blank"}
- [Single-cell best practices book](https://www.sc-best-practices.org/){:target="_blank"}
- Stark R. et al. *RNA sequencing: the teenage years.* Nature Reviews Genetics, 2019.
- Buenrostro J.D. et al. *ATAC-seq: A Method for Assaying Chromatin Accessibility Genome-Wide.* Current Protocols, 2015.
