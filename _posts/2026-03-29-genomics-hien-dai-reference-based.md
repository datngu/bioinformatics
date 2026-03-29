---
layout: post
title:  "Tư Duy Lớn Của Genomics Hiện Đại: Reference Genome Và Hệ Sinh Thái Đa Omics"
author: dat
categories: [ Bioinformatics, Genomics ]
image: assets/my_figs/ds/genomics_hien_dai.png
tags: [ genomics, theory, concept, multi-omics, overview, sequencing ]
---

Nếu phải chọn một ý tưởng duy nhất làm nền móng cho toàn bộ sinh học phân tử hiện đại, đó chính là ý tưởng về **bộ gen tham chiếu** — reference genome. Không phải một công nghệ cụ thể, không phải một phần mềm hay một công trình nghiên cứu đơn lẻ: đó là một **hệ tọa độ dùng chung** mà nhờ nó, mọi câu hỏi về ADN, RNA, protein hay cấu trúc nhiễm sắc thể đều có thể được đặt trong cùng một ngôn ngữ. Bài viết này sẽ giúp bạn hiểu tư duy này từ gốc rễ, và nhìn thấy tại sao hàng chục công nghệ khác nhau — tưởng như rời rạc — thực ra là những mảnh ghép của cùng một bức tranh.

## 1. Reference Genome — Khi Khoa Học Cần Một Ngôn Ngữ Chung

Trước khi Dự án Bộ Gen Người (Human Genome Project) hoàn thành vào năm 2003, mỗi phòng thí nghiệm nghiên cứu một đoạn ADN nhỏ một cách cô lập. Giống như hàng nghìn người mù mỗi người đang cầm một mảnh bản đồ khác nhau, nhưng không ai biết mảnh của mình ứng với vị trí nào trên tổng thể.

Reference genome thay đổi tất cả điều đó. Đây là một **bản đồ hệ gen chuẩn** — một trình tự ADN hoàn chỉnh, được lắp ráp, đánh số vị trí và công bố công khai cho từng loài. Đối với người, phiên bản hiện tại được gọi là GRCh38 (hay hg38), bao gồm khoảng 3,2 tỷ cặp base trải dài trên 23 cặp nhiễm sắc thể.

Giá trị của reference genome không nằm ở chỗ nó "hoàn hảo" hay "đại diện cho mọi người". Giá trị của nó nằm ở chỗ nó là **điểm quy chiếu dùng chung**. Khi một nhà nghiên cứu ở Hà Nội và một nhà nghiên cứu ở Boston cùng nói đến "vị trí 179,148,114 trên nhiễm sắc thể 17", họ biết chính xác mình đang nói về cùng một chỗ trên bản đồ. Đây là nền tảng của mọi hợp tác và so sánh trong sinh học phân tử hiện đại.

![Hệ tọa độ reference genome và các lớp omics]({{ site.baseurl }}/assets/my_figs/ds/genomics_hien_dai.png)
*Reference genome đóng vai trò là hệ trục tọa độ chung, nơi mọi lớp thông tin omics đều được neo vào.*

---

## 2. Tại Sao Cần Nhiều Công Nghệ? — Câu Hỏi Quyết Định Công Nghệ

Bản thân reference genome chỉ là một trình tự ADN tĩnh — nó không tự nói lên điều gì về một cá thể cụ thể, một tế bào ung thư, hay một con người đang bệnh. Để rút ra tri thức, ta cần **đo lường** và **so sánh** với nó, và mỗi loại câu hỏi đòi hỏi một loại phép đo khác nhau.

Đây chính là lý do genomics hiện đại không phải là một công nghệ mà là một **hệ sinh thái công nghệ**. Hãy nghĩ về nó như việc nghiên cứu một thành phố: bản đồ đường phố (reference genome) là cơ sở, nhưng để hiểu cuộc sống trong thành phố đó, bạn cần dữ liệu giao thông (lưu lượng RNA), dữ liệu xây dựng (cấu trúc chromatin), dữ liệu dân số (loại tế bào), và dữ liệu sản xuất kinh tế (protein và chất chuyển hóa).

Mỗi loại câu hỏi sinh học dẫn đến một lớp thông tin khác nhau:

| Câu hỏi | Lớp thông tin | Tên ngành |
|---------|---------------|-----------|
| ADN của cá thể này khác gì với tham chiếu? | Biến thể di truyền | **Genomics** |
| Gen nào đang được biểu hiện, ở mức độ nào? | Thông tin RNA | **Transcriptomics** |
| Vùng ADN nào đang "mở" hay "đóng" trong tế bào? | Trạng thái chromatin | **Epigenomics** |
| Protein nào thực sự hiện diện và ở nồng độ bao nhiêu? | Thông tin protein | **Proteomics** |
| Các phân tử nhỏ nào đang lưu hành trong tế bào? | Chất chuyển hóa | **Metabolomics** |

---

## 3. Hai Thế Hệ Đọc Bộ Gen — Và Triết Lý Đằng Sau

Để "đo lường" bất kỳ lớp thông tin nào ở trên, ta đều phải giải câu đố tương tự: làm thế nào đọc được các phân tử sinh học một cách tin cậy và có hệ thống? Câu trả lời phổ biến nhất là **giải trình tự ADN** — chuyển đổi phân tử thành dữ liệu số.

### Short-read sequencing — Sức mạnh của quy mô

Công nghệ Illumina, đại diện tiêu biểu của thế hệ này, đọc các đoạn ADN ngắn (100–300 nucleotide mỗi đoạn) nhưng với số lượng cực kỳ lớn — hàng trăm triệu đoạn trong một lần chạy. Triết lý ở đây là **độ sâu bù cho độ dài**: từng đoạn riêng lẻ ngắn, nhưng khi chồng lấp hàng triệu đoạn lên nhau, ta có thể xác định trình tự chính xác với độ tin cậy rất cao.

Điều này tạo ra một điểm mạnh và một điểm yếu rõ rệt. Điểm mạnh: độ chính xác gần như tuyệt đối (~99.9%), chi phí thấp, thông lượng cao. Điểm yếu: đoạn ngắn không thể vượt qua được những vùng ADN lặp đi lặp lại dài (như các vùng centromere), giống như bạn muốn lắp ghép một bức tranh ghép hình nhưng tất cả các mảnh đều trông giống nhau.

### Long-read sequencing — Đọc toàn bộ câu chuyện

PacBio và Oxford Nanopore đại diện cho triết lý ngược lại: thay vì đọc nhiều đoạn ngắn, hãy đọc ít đoạn hơn nhưng **mỗi đoạn dài hơn nhiều** — có thể lên đến hàng chục nghìn, thậm chí hàng triệu nucleotide.

Hãy nghĩ đến sự khác biệt giữa đọc tin nhắn và đọc tiểu thuyết: short-read giống như nhận một nghìn mảnh vụn của nhiều trang sách khác nhau; long-read giống như đọc từng chương liên tục. Đối với những vùng bộ gen phức tạp, lặp lại, hoặc khi cần hiểu cấu trúc tổng thể, long-read là không thể thay thế.

Một khả năng nữa của long-read, đặc biệt là Oxford Nanopore, là **đọc trực tiếp các tín hiệu hóa học** trên ADN mà không cần biến đổi trước — điều này cho phép phát hiện methylation (một dạng biến đổi epigenetic) ngay trong quá trình giải trình tự, thay vì phải làm thí nghiệm riêng.

Trong thực hành, hai thế hệ này **bổ sung cho nhau** nhiều hơn là cạnh tranh: long-read cung cấp khung xương cấu trúc, short-read cung cấp độ chính xác chi tiết.

---

## 4. Các Lớp Omics — Cùng Một Tọa Độ, Nhiều Chiều Thông Tin

Đây là điểm cốt lõi của toàn bộ bài viết này. Mỗi lớp omics không phải là một thế giới riêng biệt — chúng đều **chia sẻ cùng một hệ tọa độ**: vị trí trên reference genome. Điều này biến chúng từ những mảnh thông tin rời rạc thành những lớp của cùng một thực thể.

### Genomics — Bản đồ biến thể di truyền

Mỗi người trong chúng ta khác với reference genome khoảng 4–5 triệu vị trí (trong tổng số 3,2 tỷ). Những khác biệt này gọi là **biến thể di truyền**. Một số biến thể lành tính — chỉ là sự đa dạng bình thường giữa người với người. Một số khác, đặc biệt khi xuất hiện ở tế bào ung thư mà không có trong tế bào lành, là tín hiệu quan trọng về cơ chế bệnh.

Việc so sánh bộ gen của tế bào ung thư với tế bào lành từ cùng một người cho phép nhận diện các **đột biến somatic** — những biến đổi ADN chỉ xuất hiện trong khối u, không di truyền, nhưng là động cơ thúc đẩy tế bào phân chia không kiểm soát.

### Transcriptomics — Lắng nghe bộ gen đang "nói gì"

ADN là bản thảo, RNA là thông điệp đang được đọc to. Không phải mọi gen đều hoạt động ở mọi lúc và mọi tế bào — bộ gen có tính **ngữ cảnh**. Tế bào cơ tim biểu hiện một tập hợp gen khác tế bào thần kinh, dù cả hai có cùng bộ ADN.

RNA-seq đo lường **lượng RNA** được tạo ra từ mỗi gen tại một thời điểm cụ thể — đây là ảnh chụp nhanh về trạng thái "biểu hiện" của toàn bộ bộ gen. Khi so sánh tế bào bệnh và tế bào lành, ta có thể thấy gen nào đang bị "bật" quá mức, gen nào đang "tắt" bất thường, và từ đó suy ra con đường sinh học nào đang bị rối loạn.

**Single-cell RNA-seq (scRNA-seq)** đẩy khái niệm này lên một tầm cao mới: thay vì lấy trung bình tín hiệu của hàng triệu tế bào trong một mẫu mô, ta đọc biểu hiện gen của **từng tế bào riêng lẻ**. Điều này lần đầu tiên cho phép chúng ta nhìn thấy sự đa dạng ẩn bên trong một mô — phát hiện ra những tiểu quần thể tế bào hiếm, theo dõi quá trình biệt hóa, hoặc xác định tế bào nào trong khối u đang kháng thuốc.

### Epigenomics — Tầng điều khiển phía trên ADN

ADN không nằm trần trụi trong tế bào. Nó cuộn chặt quanh các protein gọi là **histone**, tạo thành cấu trúc bao gói gọi là chromatin. Trạng thái của lớp bao gói này — mở hay đóng, được đánh dấu hóa học như thế nào — quyết định gen nào có thể được đọc và gen nào bị "khóa" lại.

Đây chính là tầng **epigenomics**: những biến đổi không thay đổi trình tự ADN nhưng thay đổi cách bộ gen được sử dụng. Hai tế bào với cùng một bộ ADN có thể hoạt động hoàn toàn khác nhau vì epigenome của chúng khác nhau.

Các công nghệ như **ATAC-seq** xác định vùng nào của chromatin đang mở (tức là đang hoạt động điều hòa), **ChIP-seq** tìm nơi các protein điều hòa bắt vào ADN, còn **WGBS** đo mức độ methyl hóa (một loại "nhãn dán" hóa học) trên từng cytosine trong toàn bộ bộ gen. Tất cả những tín hiệu này đều được ánh xạ về cùng một tọa độ trên reference genome.

Đặc biệt quan trọng: **epigenome có thể thay đổi theo môi trường và bệnh lý**, trong khi genome thì không (ngoại trừ đột biến). Điều này khiến epigenomics trở thành cầu nối giữa di truyền học và môi trường sống.

### Hi-C và cấu trúc 3D — Khi bộ gen có chiều sâu

Bộ gen không phải là một sợi chỉ thẳng. Trong không gian tế bào, ADN gấp lại thành những cấu trúc ba chiều phức tạp, và vị trí không gian này có vai trò quan trọng trong điều hòa gen. Hai vùng ADN cách xa nhau hàng triệu base pair trên trình tự có thể đứng sát nhau trong không gian và tương tác với nhau để bật/tắt gen.

Công nghệ **Hi-C** chụp ảnh những tương tác này — một bản đồ về ai đang "gặp gỡ" ai trong không gian tế bào nhân. Kết quả là chúng ta hiểu tại sao một đột biến ở vùng này lại ảnh hưởng đến gen ở vùng khác tưởng như "xa" trên trình tự thẳng.

---

## 5. Multi-Omics — Khi Tất Cả Lớp Cùng Kể Một Câu Chuyện

Sức mạnh thực sự của genomics hiện đại không đến từ từng lớp omics riêng lẻ, mà từ khả năng **tích hợp nhiều lớp** để trả lời những câu hỏi mà không một lớp nào có thể trả lời một mình.

Hãy lấy ví dụ về ung thư vú: một khối u kháng thuốc điều trị. Với chỉ genomics, ta thấy có đột biến trong gen PIK3CA — nhưng đột biến này tồn tại trong cả tế bào nhạy thuốc và kháng thuốc. Transcriptomics cho thấy con đường PI3K/AKT bị kích hoạt mạnh hơn trong nhóm kháng thuốc — nhưng điều gì đang kích hoạt nó? ATAC-seq tiết lộ một vùng enhancer gần gen kháng thuốc đang "mở" ra trong tế bào kháng thuốc — nhưng tại sao? ChIP-seq xác nhận một yếu tố phiên mã đang bám vào vùng đó. Và scRNA-seq phát hiện rằng chỉ khoảng 5% tế bào trong khối u — một tiểu quần thể ẩn — mang toàn bộ đặc điểm kháng thuốc này.

Không có lớp thông tin nào trong số trên, nếu đứng một mình, có thể vẽ ra bức tranh đầy đủ. Chính nhờ chúng được neo vào cùng một hệ tọa độ reference genome, ta có thể chồng lấp và đối chiếu chúng để tìm ra cơ chế.

---

## 6. Tư Duy Từ Câu Hỏi Sinh Học — Không Phải Từ Công Nghệ

Một điều quan trọng cần tránh khi mới bắt đầu với tin sinh học là **tư duy từ công nghệ**: "Tôi có dữ liệu RNA-seq, tôi nên làm gì với nó?" Tư duy đúng đắn hơn là **từ câu hỏi sinh học**: "Tôi muốn hiểu tại sao những tế bào này hoạt động khác nhau — dữ liệu nào sẽ giúp tôi trả lời điều đó?"

Câu hỏi sinh học quyết định lớp thông tin cần thiết. Lớp thông tin quyết định công nghệ. Công nghệ quyết định quy trình phân tích. Đây là chuỗi tư duy nên đi từ trên xuống, không phải từ dưới lên.

Reference genome là nền tảng bất biến trong chuỗi này. Dù bạn đang nghiên cứu bộ gen của người, chuột, lúa hay vi khuẩn đường ruột — nguyên lý đều như nhau: xây dựng bản đồ tham chiếu, rồi đo lường mọi thứ **tương đối so với bản đồ đó**.

---

## Kết Luận

Genomics hiện đại không phải là một tập hợp công nghệ ngẫu nhiên. Đó là một **kiến trúc tư duy** được xây dựng xoay quanh một trung tâm duy nhất: reference genome. Mỗi công nghệ omics — từ giải trình tự DNA đến đo trạng thái chromatin, từ đọc RNA đến vẽ bản đồ cấu trúc 3D — là một cách nhìn vào cùng một thực thể từ một góc độ khác nhau. Và nhờ tất cả đều dùng chung một hệ tọa độ, ta có thể tích hợp chúng lại để hiểu những câu hỏi phức tạp nhất của sinh học — từ cơ chế ung thư đến tiến hóa, từ phát triển phôi đến bệnh hiếm gặp.

Trong các bài viết tiếp theo, chúng ta sẽ đi sâu vào từng lớp omics riêng biệt: cơ chế hoạt động của RNA-seq, logic đằng sau ATAC-seq, và triết lý của single-cell genomics. Những bài viết hướng dẫn thực hành (tutorial) sẽ đồng hành để giúp bạn tự tay phân tích dữ liệu thực tế.

**Tài liệu tham khảo:**
- Claussnitzer M. et al. *A brief history of human disease genetics.* Nature, 2020. — Tổng quan xuất sắc về hành trình từ gen đến bệnh.
- Stark R. et al. *RNA sequencing: the teenage years.* Nature Reviews Genetics, 2019.
- Buenrostro J.D. et al. *ATAC-seq: A Method for Assaying Chromatin Accessibility Genome-Wide.* Current Protocols, 2015.
- Lieberman-Aiden E. et al. *Comprehensive mapping of long-range interactions reveals folding principles of the human genome.* Science, 2009. — Bài gốc về Hi-C.
- [ENCODE Project](https://www.encodeproject.org/){:target="_blank"} — Kho dữ liệu epigenomics lớn nhất thế giới.
- [nf-core pipelines](https://nf-co.re/pipelines){:target="_blank"} — Pipeline phân tích chuẩn hoá cho mọi lớp omics.

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
