---
layout: post
title:  "Từ Dữ Liệu Thô Đến Biến Thể Di Truyền: Genome Alignment, Variant Calling, QC và Annotation"
author: dat
categories: [ Bioinformatics, Genomics ]
image: assets/my_figs/ds/genome-alignment-variant-calling-cover.png
tags: [ theory, concept, overview, genomics, variant-calling, alignment, annotation ]
---

Khi một mẫu DNA được đưa vào máy giải trình tự, kết quả trả về không phải là một câu trả lời hoàn chỉnh về bộ genome của cá thể đó. Nó là hàng trăm triệu đoạn ngắn rời rạc — những mảnh vỡ nhỏ của một câu chuyện lớn hơn nhiều. Việc khâu ghép những mảnh vỡ đó lại, tìm ra chỗ nào khác với tham chiếu chuẩn, loại bỏ những khác biệt giả tạo do nhiễu kỹ thuật, rồi giải thích ý nghĩa sinh học của từng biến thể — đó là bốn chương của quy trình phân tích genomics hiện đại. Bài viết này lý giải trực giác đằng sau từng chương đó, không phải từ góc độ lệnh chạy hay tham số phần mềm, mà từ góc độ tại sao chúng ta làm vậy và vấn đề khó ở đâu.

## 1. Alignment

### 1.1. Bài Toán Cốt Lõi

**Read alignment** (căn chỉnh dữ liệu giải trình tự) là bước đặt từng đoạn ngắn — gọi là *read* — vào đúng vị trí của nó trên reference genome. Hãy hình dung bạn xé nát một quyển sách thành hàng trăm triệu mảnh giấy nhỏ, mỗi mảnh chứa vài chục ký tự. Bài toán alignment là: cho mỗi mảnh giấy đó, hãy tìm xem nó xuất phát từ trang nào, dòng nào của quyển sách. Nghe đơn giản, nhưng ở quy mô genome người — 3,2 tỷ ký tự — với hơn 100 triệu mảnh nhỏ, đây là một trong những bài toán tính toán phức tạp nhất trong tin sinh học.

Khó khăn đầu tiên là **quy mô**. Nếu bạn thử tìm từng mảnh một cách ngây thơ bằng cách đọc xuyên suốt quyển sách từ đầu đến cuối, thời gian tính toán sẽ là không tưởng. Vì vậy, mọi phần mềm alignment hiện đại đều xây dựng trước một **chỉ mục** (index) của reference genome — tương tự như mục lục của một quyển từ điển hoặc chỉ mục của một cơ sở dữ liệu — để có thể "tra nhanh" vị trí gần đúng, rồi mới khớp chính xác.

### 1.2. Hai Chiến Lược Lập Chỉ Mục Quan Trọng

Hai phương pháp lập chỉ mục chiếm ưu thế trong lịch sử phát triển của các công cụ alignment là **hash-based** và **BWT-based**.

Chiến lược **hash-based** (dựa trên bảng băm) chia reference genome thành các từ khóa ngắn (*k-mers*) và lưu trữ vị trí của chúng trong một bảng băm. Khi cần tìm một read mới, ta trích xuất k-mers từ read đó rồi tra bảng — thao tác này gần như tức thì. Công cụ **ELAND** và các phiên bản đầu của **Novoalign** dùng phương pháp này. Nhược điểm là bảng băm tiêu thụ rất nhiều bộ nhớ khi k-mer kích thước lớn.

Chiến lược **BWT-based** (dựa trên Burrows-Wheeler Transform) tinh tế hơn. BWT là một phép biến đổi lý thuyết thông tin ban đầu được thiết kế để nén dữ liệu, nhưng các nhà tin sinh học nhận ra rằng nó có thể được đảo ngược sao cho việc tìm kiếm chuỗi trở nên cực kỳ hiệu quả cả về tốc độ lẫn bộ nhớ. Các công cụ như **BWA-MEM** và **Bowtie2** đều dựa trên nền tảng này. Kết hợp với cấu trúc dữ liệu **FM-index**, BWA-MEM trở thành công cụ chuẩn cho hầu hết các quy trình WGS (whole genome sequencing) và WES (whole exome sequencing) hiện nay.

### 1.3. Hơn Trăm Triệu Quyết Định Căn Chỉnh

Mỗi read đặt ra một quyết định: nó thuộc về đâu trên reference genome? Quyết định này không phải lúc nào cũng rõ ràng. Genome người chứa nhiều **vùng lặp lại** (*repetitive regions*) — những đoạn trình tự giống nhau xuất hiện ở nhiều nơi. Khi một read rơi vào vùng lặp, có thể có hàng chục vị trí khớp ngang nhau, và phần mềm phải quyết định xử lý thế nào.

Đây là lý do alignment không chỉ trả về vị trí mà còn trả về một **mapping quality score** (điểm chất lượng căn chỉnh) — một con số cho biết mức độ tin cậy của quyết định đặt read vào vị trí đó. Một read khớp duy nhất và hoàn hảo sẽ có mapping quality cao; một read rơi vào vùng lặp và có thể khớp ở nhiều chỗ sẽ có mapping quality thấp. Điểm số này sẽ được sử dụng lại trong các bước sau như một bộ lọc chất lượng.

Ngoài ra, alignment còn phải xử lý **indels** (insertion và deletion — những vị trí mà genome của cá thể có một đoạn dư ra hoặc thiếu so với reference) bằng kỹ thuật **gapped alignment** (căn chỉnh có khoảng trống). Điều này khó hơn nhiều so với khớp đơn thuần vì phần mềm phải vừa tìm vị trí vừa suy luận ra cấu trúc của indel tiềm năng.

Kết quả của alignment là một file **BAM** — một định dạng nhị phân nén chứa toàn bộ thông tin căn chỉnh của mọi read, sẵn sàng cho bước tiếp theo.

---

## 2. Variant Calling

### 2.1. Từ Sự Khác Biệt Tín Hiệu Đến Biến Thể

Sau khi alignment hoàn tất, mỗi vị trí trên reference genome được phủ bởi nhiều read chồng lên nhau. **Variant calling** (gọi biến thể di truyền) là quá trình quan sát tất cả các read tại mỗi vị trí và hỏi: *liệu sự khác biệt tôi thấy ở đây có phải là biến thể thật hay chỉ là lỗi kỹ thuật?*

Đây là một bài toán phân biệt tín hiệu khỏi nhiễu theo nghĩa đen. Mỗi lần máy đọc trình tự một nucleotide, có xác suất nhỏ nó đọc sai. Những lỗi ngẫu nhiên này trải đều và không có mẫu hình. Nhưng nếu tại một vị trí cụ thể, một phần lớn các read đều đồng loạt báo cáo cùng một nucleotide khác với reference, đó là tín hiệu đáng tin cậy rằng cá thể đó thực sự mang một **variant** (biến thể di truyền) tại vị trí đó.

### 2.2. Mô Hình Thống Kê Dưới Lớp Phần Mềm

Công cụ phổ biến nhất trong lĩnh vực này là **GATK** (*Genome Analysis Toolkit*), được phát triển bởi Broad Institute. Thuật toán cốt lõi của GATK cho variant calling là **HaplotypeCaller**, và điểm mấu chốt là: nó không nhìn từng vị trí một cách cô lập mà nhìn cả một vùng cục bộ xung quanh vị trí đó, tái dựng lại các **haplotype** (haplotype: chuỗi biến thể di truyền được kế thừa cùng nhau trên cùng một nhiễm sắc thể) có thể có, rồi mới gọi biến thể.

Cụ thể hơn, HaplotypeCaller hoạt động theo ba bước khái niệm: đầu tiên nó xác định vùng nào của genome có khả năng chứa biến thể (dựa trên sự khác biệt so với reference), sau đó với mỗi vùng đó nó xây dựng một **de Bruijn graph** (đồ thị de Bruijn) — một cấu trúc dữ liệu mô tả tất cả các haplotype khả thi từ dữ liệu đọc trình tự — rồi dùng mô hình **Bayes** để gán xác suất cho từng haplotype và quyết định biến thể nào thực sự có mặt.

Điểm mạnh của cách tiếp cận haplotype-based so với variant calling từng vị trí đơn giản là: nó xử lý tốt hơn các indels phức tạp và các biến thể nằm gần nhau, vì chúng được xét đồng thời trong cùng một khung haplotype thay vì độc lập với nhau.

### 2.3. Germline Và Somatic

Trong lĩnh vực genomics lâm sàng, có một phân biệt quan trọng: **germline variant** (biến thể dòng mầm) là biến thể di truyền từ cha mẹ, có trong mọi tế bào của cơ thể và có thể truyền cho thế hệ sau; còn **somatic variant** (biến thể soma) là đột biến xuất hiện trong quá trình sống ở một số tế bào nhất định — đây là cơ sở phân tử của ung thư.

Hai loại biến thể này đặt ra bài toán tính toán khác nhau hoàn toàn. Germline variant thường hiện diện ở 50% (dị hợp tử) hoặc 100% (đồng hợp tử) số reads tại một vị trí — những tỷ lệ rõ ràng, dễ phân biệt với nhiễu. Somatic variant có thể chỉ có mặt ở 5–10% reads (hoặc thậm chí thấp hơn trong khối u không đồng nhất), ranh giới giữa tín hiệu và nhiễu trở nên mờ nhạt hơn nhiều. Vì vậy, phần mềm như **Mutect2** (cũng trong GATK) được thiết kế riêng cho somatic calling với mô hình thống kê phức tạp hơn, thường yêu cầu thêm mẫu đối chứng từ mô lành của cùng cá thể để so sánh.

---

## 3. Post-Variant QC

### 3.1. Không Phải Mọi Biến Thể Được Gọi Đều Là Biến Thể Thật

Ngay cả khi variant caller hoạt động tốt nhất, output của nó vẫn chứa một tỷ lệ nhất định **false positive** — những "biến thể" thực ra là lỗi kỹ thuật hoặc artifact của quy trình thí nghiệm. Post-variant QC (kiểm soát chất lượng sau khi gọi biến thể) là bước lọc bỏ những false positive này trước khi đi vào phân tích xuôi.

Có một nghịch lý thú vị ở đây: nếu lọc quá chặt, ta loại bỏ cả biến thể thật (*false negatives*); nếu lọc quá lỏng, ta giữ lại nhiều artifact (*false positives*). Công việc QC là tìm điểm cân bằng phù hợp với mục đích nghiên cứu, và không có một ngưỡng duy nhất nào đúng cho mọi trường hợp.

### 3.2. Các Chiều Thông Tin Dùng Để Lọc

Mỗi biến thể được GATK (và các tool tương tự) báo cáo kèm theo nhiều **annotation** kỹ thuật — không phải annotation sinh học, mà là các số đo chất lượng kỹ thuật mô tả độ tin cậy của cuộc gọi đó. Một số chiều quan trọng:

- **Độ phủ** (*sequencing depth*): vị trí được phủ bởi bao nhiêu reads? Biến thể gọi từ rất ít reads (ví dụ dưới 8–10 reads) vốn dĩ không đáng tin cậy vì dữ liệu thống kê quá ít.
- **Chất lượng base**: các reads báo cáo biến thể đó có base quality cao không, hay chính xác tại vị trí đó chất lượng đọc bị thấp?
- **Mapping quality**: reads phủ vị trí đó có căn chỉnh tốt không, hay chúng đến từ những reads có mapping quality thấp (vùng lặp lại)?
- **Strand bias**: biến thể có xuất hiện đều ở cả hai chiều đọc (forward và reverse strand) không? Artifact kỹ thuật thường chỉ xuất hiện trên một chiều.
- **Allele balance**: tỷ lệ reads mang biến thể so với reads mang allele reference có phù hợp với kỳ vọng (50% cho dị hợp tử germline) không?

### 3.3. VQSR Và Hard Filtering

GATK cung cấp hai chiến lược lọc chính. **VQSR** (*Variant Quality Score Recalibration* — tái hiệu chỉnh điểm chất lượng biến thể) huấn luyện một mô hình thống kê trên một tập biến thể đã được xác nhận là đúng (*truth set*) như dbSNP hay HapMap, rồi dùng mô hình đó để đánh giá lại chất lượng của từng biến thể mới. Đây là phương pháp học máy ứng dụng, và nó hoạt động rất tốt khi có đủ dữ liệu — thường cần ít nhất vài chục mẫu để mô hình VQSR hội tụ.

Với dataset nhỏ hơn, người ta dùng **hard filtering** — đặt ngưỡng cứng cho từng annotation: ví dụ loại bỏ mọi biến thể có depth dưới 10, có mapping quality dưới 40, hoặc có strand bias quá lớn. Phương pháp này đơn giản hơn nhưng kém linh hoạt hơn VQSR vì nó không học được sự tương tác giữa các annotation.

### 3.4. Bước Tiền Xử Lý Làm Nền Tảng

Chất lượng của bước QC phụ thuộc không chỉ vào bản thân các ngưỡng lọc mà còn vào những bước tiền xử lý ở giai đoạn alignment. Hai bước đặc biệt quan trọng: **duplicate marking** (đánh dấu reads trùng lặp) và **base quality score recalibration** (BQSR).

Duplicate marking xử lý vấn đề PCR duplicates: trong quá trình chuẩn bị thư viện, một phân tử DNA có thể được khuếch đại nhiều lần bởi PCR, tạo ra nhiều reads giống hệt nhau không phải vì thực sự có nhiều phân tử riêng biệt. Nếu không loại bỏ chúng, biến thể ở những phân tử đó sẽ bị tính nhiều lần, làm méo tỷ lệ allele và thống kê downstream.

BQSR thừa nhận rằng điểm chất lượng do máy giải trình tự gán ban đầu thường không chính xác hoàn toàn — chúng có xu hướng bị ước tính sai một cách có hệ thống. BQSR dùng các vị trí known variant để học lại hàm hiệu chỉnh, sao cho base quality score cuối cùng phản ánh đúng xác suất lỗi thực tế của từng base.

---

## 4. Variant Annotation

### 4.1. Tên Và Địa Chỉ Chưa Đủ

Sau QC, ta có một danh sách biến thể sạch — mỗi biến thể là một tọa độ trên genome (nhiễm sắc thể, vị trí, allele reference, allele thay thế). Nhưng tọa độ thuần túy không mang ý nghĩa sinh học. **Variant annotation** (chú thích biến thể) là quá trình gắn thêm ngữ cảnh sinh học vào từng biến thể: nó nằm trong gene nào, ảnh hưởng đến protein thế nào, đã từng được ghi nhận chưa, và có liên quan đến bệnh gì không.

Annotation là cầu nối từ tín hiệu genomics sang tri thức sinh học y học, và nó thường đặt ra nhiều câu hỏi mới hơn là giải quyết câu hỏi cũ.

### 4.2. Phân Loại Chức Năng

Bước đầu tiên của annotation là xác định **functional consequence** (hậu quả chức năng) của biến thể dựa trên vị trí của nó trong cấu trúc gene. Công cụ như **VEP** (*Variant Effect Predictor* của Ensembl) hay **SnpEff** tự động thực hiện việc này.

Một biến thể trong vùng mã hóa (*coding region*) có thể là:
- **Synonymous variant** (biến thể đồng nghĩa): thay đổi nucleotide nhưng không thay đổi amino acid do tính thoái hóa của bảng mã di truyền — thường được xem là ít tác động nhất.
- **Missense variant** (biến thể sai nghĩa): thay đổi một amino acid, có thể vô hại hoặc có thể phá vỡ chức năng protein tùy thuộc vào amino acid đó có bảo tồn tiến hóa hay không và nằm ở vị trí cấu trúc nào.
- **Nonsense variant** (biến thể vô nghĩa): tạo ra codon dừng sớm, thường cắt ngắn và bất hoạt protein.
- **Splice-site variant** (biến thể vị trí nối exon-intron): làm rối loạn quá trình splicing, có thể dẫn đến mất exon hoặc giữ lại intron trong mRNA trưởng thành.

Biến thể ngoài vùng mã hóa — trong **promoter**, **enhancer**, **UTR** (untranslated region), hay **intron** — cũng có thể có tác động sinh học quan trọng bằng cách thay đổi gene expression regulation, nhưng chúng khó tiên đoán hơn nhiều bằng tính toán đơn thuần.

### 4.3. Cơ Sở Dữ Liệu Như Bộ Nhớ Tập Thể

Annotation mạnh mẽ không chỉ dựa vào suy luận *ab initio* (từ bước đầu tiên, không dựa trên kiến thức trước đó) mà còn tận dụng **cơ sở dữ liệu thực nghiệm** tích lũy từ hàng chục năm nghiên cứu. Một số cơ sở dữ liệu nền tảng:

| Cơ sở dữ liệu | Nội dung |
|----------------|----------|
| **dbSNP** | Biến thể đã được đánh mã rsID, thường gặp trong quần thể |
| **ClinVar** | Biến thể có bằng chứng lâm sàng liên quan đến bệnh |
| **gnomAD** | Tần số allele trong >140.000 genome và exome quần thể toàn cầu |
| **COSMIC** | Somatic variant trong ung thư từ hàng triệu khối u |
| **OMIM** | Gene–bệnh associations từ y văn |

*Bảng 4.1. Các cơ sở dữ liệu cốt lõi trong variant annotation.*

Một biến thể đã có trong dbSNP với tần số cao (ví dụ >1% trong quần thể) gần như chắc chắn là benign (*lành tính*). Một biến thể xuất hiện trong ClinVar với phân loại *pathogenic* và bằng chứng mạnh từ nhiều nguồn độc lập là căn cứ lâm sàng quan trọng. Đây là cách annotation kết hợp dữ liệu genomics cá nhân với bộ nhớ tập thể của cộng đồng khoa học.

### 4.4. Tiên Đoán Tác Động Bằng Học Máy

Không phải mọi biến thể đều đã có trong cơ sở dữ liệu — đặc biệt các biến thể hiếm hoặc mới. Với những trường hợp này, cần các công cụ **in silico prediction** (tiên đoán tính toán) để ước tính tác động chức năng từ các đặc điểm của biến thể.

**SIFT** dựa trên nguyên lý bảo tồn tiến hóa: nếu một amino acid được bảo tồn qua hàng chục triệu năm tiến hóa ở nhiều loài khác nhau, thay đổi nó nhiều khả năng sẽ gây hại cho chức năng protein. **PolyPhen-2** kết hợp bảo tồn tiến hóa với thông tin cấu trúc 3D của protein. Các công cụ thế hệ mới như **CADD** (*Combined Annotation Dependent Depletion*) tổng hợp hàng chục loại annotation thành một điểm tổng hợp duy nhất bằng phương pháp học máy, trong khi **AlphaMissense** — được xây dựng dựa trên kiến trúc của AlphaFold — có thể tiên đoán tác động của missense variant với độ chính xác gần với đánh giá thực nghiệm.

Tuy vậy, tất cả các công cụ này đều chỉ là dự đoán và có giới hạn. Chúng hữu ích để ưu tiên các biến thể cần nghiên cứu sâu hơn, không phải để kết luận cuối cùng về tính bệnh lý của một biến thể cụ thể.

### 4.5. ACMG và Phân Loại Lâm Sàng

Trong bối cảnh y học, việc annotation phải dẫn đến một **phân loại lâm sàng** để có thể đưa ra quyết định chẩn đoán. Hiệp hội Di truyền Y học Hoa Kỳ (ACMG) đề xuất một framework chuẩn hóa với năm mức phân loại: *Pathogenic* (gây bệnh), *Likely Pathogenic* (có khả năng gây bệnh), *Variant of Uncertain Significance* — VUS (biến thể không xác định ý nghĩa), *Likely Benign* (có khả năng lành tính) và *Benign* (lành tính).

Mỗi biến thể được chấm điểm dựa trên một bộ tiêu chí bằng chứng chuẩn hóa — gọi là **ACMG criteria** — kết hợp tần số quần thể, bằng chứng chức năng, dữ liệu phả hệ, và tiên đoán tính toán. Một phần đáng kể biến thể hiếm phát hiện trong lâm sàng rơi vào hạng mục VUS, tạo ra thách thức thực tế: làm thế nào để truyền đạt sự không chắc chắn này đến bệnh nhân và bác sĩ?

---

## 5. Bức Tranh Tổng Thể

Bốn giai đoạn — alignment, variant calling, QC, và annotation — không phải là bốn bước rời rạc mà là bốn lớp quyết định lồng vào nhau. Chất lượng của mỗi lớp phụ thuộc vào lớp trước, và mỗi lớp có những điểm không chắc chắn riêng của nó.

Alignment quyết định mỗi mảnh dữ liệu được "đặt" ở đâu trên bản đồ genome. Variant calling quyết định những khác biệt nào giữa mảnh đó và bản đồ là tín hiệu thật. QC loại bỏ những tín hiệu mà bằng chứng kỹ thuật không đủ mạnh. Annotation giải thích những tín hiệu còn lại trong ngôn ngữ sinh học và y học.

Điều đáng chú ý là quy trình này liên tục được cải thiện cả về hai phương diện: phần cứng (công nghệ giải trình tự thế hệ mới, đặc biệt long-read sequencing từ PacBio và Oxford Nanopore, cho phép đọc trực tiếp các vùng lặp lại và structural variant phức tạp mà short-read vẫn còn hạn chế) và phần mềm (các mô hình deep learning đang dần thay thế các phương pháp thống kê cổ điển trong alignment và variant calling).

Chúng tôi sẽ đề cập quy trình thực hành — bao gồm cài đặt công cụ, cú pháp lệnh và cách giải thích kết quả — trong bài tutorial riêng.

## Kết Luận

Từ hàng trăm triệu mảnh ngắn đến một danh sách biến thể có ý nghĩa lâm sàng là một hành trình tính toán phức tạp, nhưng logic của từng bước đều có thể hiểu được từ góc độ trực giác. Alignment giải quyết bài toán định vị quy mô lớn bằng chỉ số thông minh. Variant calling phân biệt tín hiệu khỏi nhiễu bằng mô hình haplotype-based và thống kê Bayes. Post-variant QC loại bỏ những tín hiệu thiếu bằng chứng bằng cách nhìn vào nhiều chiều kỹ thuật đồng thời. Annotation kết nối dữ liệu genomics cá nhân với kho tàng tri thức tập thể từ hàng triệu genome đã nghiên cứu trước.

Hiểu rõ trực giác đằng sau những lớp này là điều kiện tiên quyết để sử dụng các công cụ hiệu quả, diễn giải kết quả đúng đắn, và nhận ra giới hạn của phân tích — những kỹ năng không thể có được chỉ từ việc học cú pháp của một phần mềm cụ thể.

## Tài Liệu Tham Khảo

Li, H., & Durbin, R. (2009). Fast and accurate short read alignment with Burrows-Wheeler Aligner.
*Bioinformatics*, *25*(14), 1754–1760. https://doi.org/10.1093/bioinformatics/btp324

McKenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky, A., Garimella, K., Altshuler, D., Gabriel, S., Daly, M., & DePristo, M. A. (2010). The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data.
*Genome Research*, *20*(9), 1297–1303. https://doi.org/10.1101/gr.107524.110

Poplin, R., Ruano-Rubio, V., DePristo, M. A., Fennell, T. J., Carneiro, M. O., Van der Auwera, G. A., Kling, D. E., Gauthier, L. D., Levy-Moonshine, A., Roazen, D., Shakir, K., Thibault, J., Chandramouliswaran, I., Gabriel, S., Altshuler, D., Daly, M., & Bank, E. (2018). Scaling accurate genetic variant discovery to tens of thousands of samples.
*bioRxiv*. https://doi.org/10.1101/201178

McLaren, W., Gil, L., Hunt, S. E., Riat, H. S., Ritchie, G. R. S., Thormann, A., Flicek, P., & Cunningham, F. (2016). The Ensembl Variant Effect Predictor.
*Genome Biology*, *17*(1), 122. https://doi.org/10.1186/s13059-016-0974-4

Richards, S., Aziz, N., Bale, S., Bick, D., Das, S., Gastier-Foster, J., Grody, W. W., Hegde, M., Lyon, E., Spector, E., Voelkerding, K., & Rehm, H. L. (2015). Standards and guidelines for the interpretation of sequence variants: A joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology.
*Genetics in Medicine*, *17*(5), 405–424. https://doi.org/10.1038/gim.2015.30

Cheng, J., Novati, G., Pan, J., Bycroft, C., Žemgulytė, A., Applebaum, T., Pritzel, A., Wong, L. H., Zuk, O., Jumper, J., Hassabis, D., Kohli, P., & Žemgulytė, A. (2023). Accurate proteome-wide missense variant effect prediction with AlphaMissense.
*Science*, *381*(6664), eadg7492. https://doi.org/10.1126/science.adg7492
