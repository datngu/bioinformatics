---
layout: post
title: "Gene Regulatory Network: Từ Phân Tử Đến Multi-Omics Và Ứng Dụng Lâm Sàng"
author: dat
categories: [ Bioinformatics, Genomics ]
image: assets/my_figs/ds/gene-regulatory-network.png
tags: [ genomics, bioinformatics, theory, concept, overview, gene-regulation, multi-omics, single-cell ]
---

Genome của người chứa khoảng 20.000 gene mã hóa protein, nhưng không phải gene nào cũng hoạt động ở mọi tế bào trong mọi thời điểm. Một tế bào gan, một tế bào thần kinh, và một tế bào miễn dịch đều mang cùng một bộ DNA, nhưng chúng có hình dạng, chức năng, và tuổi thọ hoàn toàn khác nhau. Bí mật nằm ở **gene regulatory network** (GRN: mạng lưới điều hòa gene, hệ thống các phân tử và tín hiệu quy định việc gene nào được bật, gene nào bị tắt, và ở mức độ nào trong từng loại tế bào tại từng thời điểm cụ thể). GRN là lớp điều khiển cao nhất của sinh học tế bào — nó quyết định **bản sắc tế bào**, phản ứng với môi trường, và cả quá trình phát bệnh khi mạng lưới này bị rối loạn.

Hiểu GRN không chỉ là câu hỏi cơ bản của sinh học phát triển mà còn là chìa khóa để giải mã cơ chế ung thư, tìm drug target mới, và điều khiển việc tái lập trình tế bào. Bài viết này dẫn dắt từ các thành phần phân tử cơ bản nhất của GRN đến các phương pháp multi-omics hiện đại và các ứng dụng thực tiễn vào năm 2026.

---

## 1. Khái Niệm Nền Tảng

### 1.1. Định Nghĩa Và Cấu Trúc Mạng Lưới

Một GRN về bản chất là một **đồ thị có hướng** (directed graph) trong đó mỗi **node** đại diện cho một gene (hoặc sản phẩm của nó), và mỗi **edge** (cạnh) được gán hướng và dấu — dương cho kích hoạt, âm cho ức chế. Nếu **transcription factor** (TF: nhân tố phiên mã, là protein gắn vào DNA để bật hoặc tắt quá trình sản xuất RNA từ một gene) A kích hoạt biểu hiện gene B, ta vẽ một mũi tên từ A đến B. Nếu TF C ức chế gene D, mũi tên từ C đến D mang dấu âm.

Trong thực tế, GRN không đơn giản chỉ là các mối quan hệ TF–gene đơn lẻ. Các gene mã hóa TF, các TF đó lại điều hòa gene khác bao gồm cả các TF khác, tạo nên những **vòng phản hồi** (feedback loop) và **vòng truyền tiếp** (feedforward loop) phức tạp. Chính những vòng lặp này tạo nên tính ổn định và tính quyết định của bản sắc tế bào, cho phép tế bào "nhớ" trạng thái của mình ngay cả khi tín hiệu gốc không còn nữa.

### 1.2. Tại Sao GRN Quan Trọng

Khi các nhà khoa học giải mã genome người năm 2003, điều bất ngờ nhất không phải là 20.000 gene protein-coding — mà là **98% genome còn lại** không mã hóa protein. Phần lớn vùng DNA này, từng bị gọi là "junk DNA", thực ra chứa hàng triệu **cis-regulatory element** (CRE: yếu tố điều hòa cis, vùng DNA không mã hóa nhưng điều phối việc bật/tắt gene lân cận bằng cách gắn kết các TF và co-activator). Phát hiện này lật ngược quan điểm cũ: thay vì "gene là tất cả", sinh học hiện đại nhận ra rằng **cách gene được điều hòa quan trọng không kém cấu trúc của gene đó**.

Bằng chứng lâm sàng rõ ràng nhất đến từ các nghiên cứu genome-wide association study (GWAS: nghiên cứu liên kết toàn genome tìm các biến thể di truyền liên quan đến bệnh lý). Hơn 80% các biến thể di truyền liên quan đến bệnh — kể cả ung thư, tiểu đường, bệnh tim — nằm trong các vùng **không mã hóa protein**. Chúng không phá vỡ cấu trúc protein, nhưng chúng phá vỡ CRE, tức là phá vỡ một nút trong GRN. Đây là lý do tại sao giải mã GRN trở thành mục tiêu trung tâm của y học phân tử hiện đại.

---

## 2. Các Thành Phần Phân Tử

### 2.1. Transcription Factors

TF là trung tâm của mọi GRN. Protein này chứa một đầu **DNA-binding domain** (DBD: vùng gắn DNA, nhận diện một trình tự DNA ngắn đặc trưng gọi là motif, thường dài 6–12 nucleotide) và một đầu **activation/repression domain** tương tác với bộ máy phiên mã hoặc co-repressor. Genome người mã hóa khoảng **1.639 TF** (Lambert và cộng sự, 2018), mỗi TF nhận diện một tập hợp CRE khác nhau trên genome.

Điều làm cho TF thú vị về mặt mạng lưới là **tính đặc hiệu theo ngữ cảnh** (context-specificity). Cùng một TF — ví dụ **GATA1** — có thể kích hoạt bộ gene này trong tế bào hồng cầu nhưng kích hoạt bộ gene hoàn toàn khác trong tế bào khổng lồ megakaryocyte. Sự khác biệt đến từ các **co-factor** (đối tác cùng gắn DNA), trạng thái **chromatin** tại vùng đó, và cấu trúc không gian 3D của genome. Một TF không hoạt động đơn lẻ mà luôn hoạt động trong **combinatorial logic** — sự kết hợp của 2–5 TF tại một CRE mới tạo nên quyết định cuối cùng.

### 2.2. Cis-Regulatory Elements

**Promoter** (vùng khởi động: vùng DNA ngay phía trên DNA template của gene, nơi RNA polymerase và các TF cơ bản lắp ráp để bắt đầu phiên mã) là CRE quen thuộc nhất. Nhưng nó chỉ là một loại trong một hệ thống phong phú hơn nhiều:

- **Enhancer** (tăng cường tử: vùng DNA xa gene mục tiêu, đến vài trăm kilobase, nhưng có thể tiếp xúc vật lý với promoter qua vòng DNA trong không gian 3D để kích hoạt biểu hiện gene). Một gene có thể được điều tiết bởi hàng chục enhancer khác nhau, mỗi enhancer hoạt động trong một loại tế bào hoặc điều kiện khác nhau.
- **Silencer** (trình tự im lặng: CRE khi được gắn bởi repressor TF sẽ ức chế phiên mã của gene lân cận).
- **Insulator** (trình tự cách ly: CRE ngăn tín hiệu kích hoạt từ enhancer "tràn" sang gene không phải mục tiêu bằng cách thiết lập ranh giới cấu trúc trong không gian 3D).

Những phát hiện từ dự án **ENCODE** (Encyclopedia of DNA Elements: dự án lập bản đồ toàn diện các yếu tố chức năng trong genome người) đã xác định hơn 1,2 triệu CRE tiềm năng trong genome người — gấp ~60 lần số gene protein-coding. Điều này cho thấy quy mô thực sự của bài toán điều hòa gene.

### 2.3. Kiến Trúc Không Gian Ba Chiều

GRN không thể hiểu trọn vẹn nếu chỉ xem genome như một chuỗi 1D. Trong nhân tế bào, DNA được gấp cuộn theo thứ bậc:

- **Topologically associating domain** (TAD: miền liên kết theo cấu trúc, vùng genome dài 0,1–1 megabase trong đó các đoạn DNA tương tác với nhau nhiều hơn với DNA bên ngoài miền). TAD đóng vai trò như ranh giới địa lý: một enhancer nằm trong TAD của một gene sẽ khó kích hoạt gene ở TAD khác.
- **Chromatin loop** (vòng chromatin: tương tác trực tiếp giữa hai điểm xa nhau trên genome, thường được duy trì bởi protein CTCF và cohesin). Vòng này thường kết nối enhancer với promoter của gene mục tiêu.
- **A/B compartment** (phân khu A/B: ở mức độ lớn hơn, genome được chia thành vùng A — mở, phiên mã tích cực — và vùng B — đóng, gene câm lặng). Một gene chuyển từ compartment B sang A thường đồng hành với quá trình tế bào biệt hóa.

### 2.4. Motif Mạng Lưới

Không phải mọi cấu trúc cục bộ trong GRN đều xuất hiện ngẫu nhiên. Một số **network motif** (môtif mạng: mẫu kết nối lặp đi lặp lại trong mạng nhiều hơn kỳ vọng ngẫu nhiên, có chức năng điều hòa đặc trưng) xuất hiện với tần số cao và đảm nhận chức năng tính toán quan trọng:

| Motif | Chức năng điều hòa |
|-------|-------------------|
| **Feedforward loop** (FFL) | Lọc nhiễu tín hiệu; tạo trễ hay đẩy nhanh phản ứng tùy loại FFL |
| **Autoregulation dương** | Khuếch đại tín hiệu và ổn định bản sắc tế bào qua bistability |
| **Autoregulation âm** | Giảm dao động, ổn định mức biểu hiện gene về giá trị đặt trước |
| **Mutual inhibition** | Tạo quyết định nhị phân (cell fate switch), quyết định tế bào đi theo hướng phân hóa nào |

*Bảng 2.4. Các motif mạng lưới thường gặp và ý nghĩa điều hòa của chúng.*

Mô hình **toggle switch** nổi tiếng — hai TF ức chế lẫn nhau — giải thích tại sao tế bào gốc khi biệt hóa có xu hướng chọn một trong hai hướng một cách rõ ràng thay vì dừng ở trạng thái trung gian.

---

## 3. Đo Lường GRN Bằng Công Nghệ Omics

### 3.1. ChIP-seq

**ChIP-seq** (Chromatin Immunoprecipitation sequencing: kết hợp kết tủa miễn dịch chromatin với giải trình tự thế hệ mới để tạo bản đồ genome-wide về vị trí gắn kết của một TF hoặc histone mark cụ thể) là tiêu chuẩn vàng để lập bản đồ TF–DNA binding. Nguyên lý: cố định protein–DNA bằng formaldehyde, phân mảnh chromatin, rồi dùng kháng thể đặc hiệu kéo xuống phức hợp protein–DNA, giải trình tự DNA phần bị kéo xuống, và định vị các đỉnh (peaks) để biết TF đó gắn ở đâu trên genome.

ChIP-seq cho phép suy luận mạng lưới theo chiều **từ TF đến target gene**: biết TF X gắn vào CRE của gene Y kéo theo giả thuyết X điều hòa Y. Kết hợp dữ liệu ChIP-seq của hàng chục TF, ta xây dựng được một bản đồ TF-binding toàn genome. Giới hạn của ChIP-seq là cần lượng tế bào lớn (vài triệu) và mỗi thí nghiệm chỉ lập bản đồ được **một** TF hoặc histone mark tại một thời điểm.

### 3.2. ATAC-seq

**ATAC-seq** (Assay for Transposase-Accessible Chromatin with sequencing: sử dụng enzyme transposase Tn5 để cắt và gắn adapter vào những vùng chromatin hở — tức vùng không bị nucleosome che phủ — sau đó giải trình tự để lập bản đồ chromatin accessibility toàn genome) giải quyết vấn đề khác: không xác định TF cụ thể nào gắn, mà xác định vùng DNA nào đang **mở hay đóng**. Vùng chromatin hở là dấu hiệu của CRE đang hoạt động — và từ đó có thể suy luận TF nào đang gắn vào dựa trên **motif enrichment** (phân tích motif: tìm kiếm trình tự ngắn đặc trưng của các TF trong vùng mở).

ATAC-seq đặc biệt mạnh mẽ vì chỉ cần **~500 tế bào** (so với vài triệu của ChIP-seq) và có thể thực hiện ở cấp độ single-cell (scATAC-seq). Điều này mở ra khả năng lập bản đồ CRE riêng cho từng tế bào trong một mô không đồng nhất.

### 3.3. Hi-C Và Chromosome Conformation Capture

**Hi-C** (High-throughput Chromosome Conformation Capture: kỹ thuật đo lường toàn genome về các tương tác DNA–DNA trong không gian 3D của nhân tế bào, cho phép lập bản đồ TAD, compartment, và vòng chromatin) đã cách mạng hóa hiểu biết về kiến trúc 3D genome. Bằng cách bắt và giải trình tự các đầu DNA đang tiếp xúc vật lý với nhau, Hi-C tạo ra một **contact map** — ma trận đối xứng ghi lại tần suất tương tác giữa mọi cặp vùng genome.

Từ contact map, ta xác định được: TAD (từ vùng trống xung quanh đường chéo), A/B compartment (từ tương quan của tín hiệu eigenvector), và chromatin loop (các "dots" sáng trên đường chéo biểu thị cặp vùng DNA tiếp xúc trực tiếp). Bổ sung chiều thứ ba này vào GRN quan trọng vì nó giải thích tại sao enhancer xa hàng trăm kilobase vẫn điều hòa được một gene cụ thể — chúng tiếp xúc trực tiếp qua vòng chromatin.

### 3.4. RNA-seq Và Inference Từ Biểu Hiện Gene

**RNA-seq** (giải trình tự RNA) đo lường **transcriptome** — tập hợp toàn bộ phân tử RNA được tổng hợp trong tế bào tại một thời điểm — và vẫn là trụ cột của phân tích biểu hiện gene. Từ góc độ GRN, RNA-seq cho phép:

- **Suy luận mạng lưới từ dữ liệu biểu hiện** (expression-based GRN inference): nếu TF X và gene Y luôn tăng hoặc giảm cùng nhau qua nhiều điều kiện và mẫu, có thể có quan hệ điều hòa giữa chúng. Các thuật toán như **ARACNE** (dựa trên mutual information), **GENIE3** (dựa trên random forest), và **GRNBoost2** khai thác tương quan này.
- **Thêm chiều thông tin** phân tử cho các mạng được xây dựng từ ChIP-seq hay ATAC-seq bằng cách kiểm tra xem sự hiện diện của TF binding có thực sự kéo theo thay đổi biểu hiện gene hay không.

---

## 4. Single-Cell Omics Và Suy Luận GRN

### 4.1. Tại Sao Cần Single-Cell

Phần lớn dữ liệu omics cổ điển là **bulk data** — trung bình hóa trên hàng triệu tế bào cùng lúc. Điều này che giấu **tính dị chất tế bào** (cellular heterogeneity): trong một khối u hay một mô phức tạp, hàng chục loại tế bào khác nhau tồn tại cùng nhau, mỗi loại có GRN riêng. Bulk RNA-seq cho ra một tín hiệu hỗn hợp không thể phân giải được.

**scRNA-seq** (single-cell RNA sequencing: giải trình tự RNA ở cấp độ từng tế bào đơn lẻ, cung cấp hồ sơ biểu hiện gene độc lập cho mỗi tế bào trong một mẫu) giải quyết vấn đề này. Với công nghệ droplet-based như **10x Genomics Chromium**, ta có thể đo transcriptome của hàng chục nghìn tế bào đơn lẻ trong một thí nghiệm. Kết quả là một ma trận gene x tế bào — nền tảng để xây dựng GRN riêng cho từng loại tế bào hoặc trạng thái tế bào.

### 4.2. scATAC-seq Và Single-Cell Epigenomics

**scATAC-seq** (single-cell ATAC-seq: phiên bản single-cell của ATAC-seq, lập bản đồ chromatin accessibility riêng cho từng tế bào) bổ sung chiều epigenomics vào phân tích single-cell. Sự mở/đóng của CRE ở cấp độ từng tế bào phản ánh trực tiếp trạng thái điều hòa của tế bào đó — thường nhạy hơn RNA vì thay đổi chromatin thường đi trước thay đổi biểu hiện gene.

Kết hợp scRNA-seq và scATAC-seq, nhà nghiên cứu có thể liên kết **CRE mở** (inferred TF binding) với **biểu hiện gene** ở cùng một tế bào (hay loại tế bào), từ đó suy luận "TF nào, tại vùng nào, điều hòa gene nào" ở độ phân giải single-cell.

### 4.3. Multiome: Đo Đồng Thời RNA Và Chromatin

**10x Genomics Multiome ATAC + Gene Expression** là công nghệ đo **đồng thời** scRNA-seq và scATAC-seq từ **cùng một tế bào**. Điều này loại bỏ vấn đề batch effect khi kết hợp hai thí nghiệm riêng lẻ và cho phép phân tích nhân quả trực tiếp hơn: biết chromatin accessibility của tế bào này → biết gene nào đang biểu hiện từ cùng tế bào đó.

**ArchR**, **Signac**, và **EpiScanpy** là các framework phân tích phổ biến hỗ trợ dữ liệu multiome. Năm 2025–2026, multiome đã trở thành tiêu chuẩn cho các nghiên cứu cell atlas quy mô lớn như Human Cell Atlas giai đoạn 2.

### 4.4. Công Cụ Suy Luận GRN Từ Single-Cell

Từ dữ liệu single-cell, các công cụ suy luận GRN chuyên biệt đã ra đời:

- **SCENIC** (Single-Cell rEgulatory Network Inference and Clustering: kết hợp correlation-based gene co-expression modules với motif enrichment trong vùng regulatory để xác định bộ TF + regulon tương ứng trong từng tế bào). SCENIC xây dựng **regulon** — TF cùng với tập target gene của nó — và tính điểm hoạt động của mỗi regulon ở từng tế bào.
- **CellOracle** (Kamimoto và cộng sự, 2023): kết hợp scATAC-seq để xác định mạng lưới cơ sở, sau đó dùng scRNA-seq để ước tính trọng số của từng kết nối. CellOracle đặc biệt mạnh ở khả năng **in silico perturbation** — mô phỏng điều gì xảy ra khi knock out một TF cụ thể.
- **Pando**: framework multi-omics inference tích hợp scRNA-seq, scATAC-seq, và dữ liệu Hi-C.
- **DIRECT-NET** và **FigR**: tập trung vào liên kết CRE với gene target dựa trên đồng biến thiên trong dữ liệu multiome.

---

## 5. Multi-Omics Integration Cho GRN

### 5.1. Tại Sao Một Lớp Dữ Liệu Không Đủ

Mỗi công nghệ omics nhìn vào GRN từ một góc độ riêng và mỗi góc độ đều có điểm mù:

- ChIP-seq cho biết TF gắn ở đâu, nhưng **không biết** binding đó dẫn đến kích hoạt hay ức chế.
- RNA-seq cho biết gene nào đang biểu hiện, nhưng **không biết** TF nào chịu trách nhiệm.
- ATAC-seq cho biết CRE nào đang mở, nhưng **không biết** TF cụ thể nào đang sử dụng vùng mở đó.
- Hi-C cho biết hai vùng DNA ở gần nhau trong không gian, nhưng **không biết** liên kết đó có chức năng điều hòa thực sự không.

Chỉ khi tích hợp đủ nhiều lớp, ta mới có thể suy luận một GRN chính xác và có chiều sâu cơ chế.

### 5.2. Framework Tích Hợp Đa Lớp

Về mặt khái niệm, có ba chiến lược tích hợp multi-omics cho GRN:

**Early integration** (tích hợp sớm): nối các lớp omics thành một ma trận đặc trưng duy nhất trước khi suy luận mạng. Đơn giản nhưng bị ảnh hưởng bởi sự chênh lệch quy mô và phân phối giữa các lớp.

**Late integration** (tích hợp muộn): suy luận GRN độc lập từ từng lớp, sau đó hợp nhất các mạng thành phần bằng cách cộng điểm trọng số. Linh hoạt hơn nhưng mất đi các tín hiệu cross-modal.

**Joint modeling** (mô hình hóa song hành): xây dựng một mô hình thống kê hoặc deep learning học đồng thời từ tất cả các lớp. Đây là xu hướng chủ đạo từ 2023 trở đi — các phương pháp như **MOFA+** (Multi-Omics Factor Analysis), **scVI-multi**, và **MultiVI** thuộc về nhóm này.

**LINGER** (Larsson và cộng sự, 2024) là ví dụ tiêu biểu: một framework deep learning học GRN từ dữ liệu bulk multi-omics quy mô lớn, sau đó dùng single-cell ATAC + RNA để cá nhân hóa GRN xuống cấp độ tế bào. Kết quả là một GRN có độ phân giải chưa từng đạt được với các phương pháp cổ điển.

### 5.3. Spatial Multi-Omics Và Vị Trí Không Gian Mô

Năm 2025–2026, **spatial transcriptomics** (đo biểu hiện gene đồng thời giữ nguyên vị trí không gian 2D của mỗi điểm trong mô) kết hợp với **spatial ATAC-seq** đã cho phép khám phá một chiều mới: GRN không chỉ khác nhau giữa loại tế bào mà còn khác nhau tùy **vị trí của tế bào trong mô**. Tế bào ở trung tâm khối u có GRN khác với tế bào ở rìa — và sự khác biệt đó liên quan trực tiếp đến tiên lượng và kháng thuốc.

**Slide-tags**, **MERFISH**, và **Stereo-seq** là các nền tảng spatial omics được xuất bản từ 2023–2025, đang mở ra kỷ nguyên **spatial GRN atlas** — bản đồ GRN theo cả loại tế bào lẫn vị trí không gian mô.

---

## 6. Deep Learning Và Foundation Models Cho GRN

### 6.1. Graph Neural Networks

**Graph neural network** (GNN: mạng nơ-ron được thiết kế để học trên dữ liệu cấu trúc đồ thị, lan truyền thông tin qua các cạnh kết nối để học biểu diễn của từng node trong ngữ cảnh với các node láng giềng) là kiến trúc phù hợp tự nhiên cho GRN vì bản thân GRN đã là một đồ thị. Các mô hình như **GRNBoost2-GNN** và **DeepDRIM** dùng GNN để suy luận cạnh trong GRN từ dữ liệu biểu hiện gene, cho phép nắm bắt được các quan hệ gián tiếp mà correlation đơn giản không bắt được.

### 6.2. Sequence-To-Function Models

Một hướng tiếp cận khác là học trực tiếp từ **trình tự DNA thô** để dự đoán hoạt động của CRE. **Basenji** và **Enformer** (Avsec và cộng sự, 2021) là các mô hình CNN + Transformer học ánh xạ từ chuỗi DNA 200 kb → profile biểu hiện gene hay ATAC-seq signal ở nhiều loại tế bào đồng thời. Enformer đạt độ chính xác đủ cao để **dự đoán tác động của đột biến** tại CRE lên biểu hiện gene — lần đầu tiên kết nối trực tiếp biến thể di truyền với GRN.

Năm 2025, **Borzoi** cải tiến Enformer bằng cách mở rộng receptive field lên 524 kb và bổ sung dự đoán RNA splicing, cho kết quả tốt hơn đáng kể trong việc mô hình hóa các enhancer xa.

### 6.3. Foundation Models Cho Single-Cell

**Geneformer** (Theodoris và cộng sự, 2023) là large language model được huấn luyện trên 29,9 triệu tế bào người từ scRNA-seq, học cách biểu diễn tế bào và gene trong cùng một không gian embedding. Ứng dụng GRN của Geneformer là **in silico perturbation**: dự đoán hiệu quả của việc knock down một TF lên transcriptome toàn tế bào, với hiệu quả tương đương thí nghiệm thực nhưng nhanh hơn vài bậc.

**scGPT** (Cui và cộng sự, 2024) và **scFoundation** (Hao và cộng sự, 2024) mở rộng paradigm này với kiến trúc GPT/BERT phi tuyến, có khả năng tích hợp đa omics trong cùng một mô hình. Năm 2026, các foundation model này đang được thử nghiệm trong **pipeline dự đoán drug target** bằng cách kết hợp in silico perturbation với dữ liệu lâm sàng.

### 6.4. Hạn Chế Của Các Mô Hình Hiện Tại

Dù ấn tượng, các deep learning model cho GRN vẫn đối mặt với các thách thức cơ bản. Đầu tiên là **interpretability** (khả năng giải thích): một mô hình dự đoán tốt nhưng không cho biết cơ chế nào đang hoạt động có giá trị hạn chế trong nghiên cứu cơ bản. Thứ hai là sự phụ thuộc vào **dữ liệu huấn luyện** — nếu tế bào hay bệnh lý cần dự đoán khác xa phân phối huấn luyện, mô hình có thể "hallucinate" quan hệ không tồn tại. Kết hợp giữa **mechanistic model** (mô hình cơ chế) và deep learning đang là xu hướng nghiên cứu tích cực để giải quyết vấn đề này.

---

## 7. Ứng Dụng Lâm Sàng Và Nghiên Cứu

### 7.1. GRN Trong Ung Thư

Ung thư là bệnh lý GRN ở cấp độ cơ bản nhất. Các **master TF** (TF chủ chốt điều phối chương trình phiên mã toàn diện của một loại tế bào) như MYC, TP53, RUNX1, và SOX2 bị đột biến, khuếch đại, hay mất chức năng trong hầu hết các loại ung thư. Khi một master TF thay đổi, toàn bộ GRN mà nó điều phối bị đảo lộn — đây là lý do tại sao một đột biến điểm duy nhất có thể dẫn đến sự thay đổi biểu hiện của hàng trăm gene.

**Super-enhancer** (vùng cis-regulatory đặc biệt mạnh, là sự kết tụ của nhiều enhancer thông thường tại một locus, gắn kết mật độ cao các TF và co-activator) đặc biệt liên quan đến ung thư. Tế bào ung thư thường tái lập trình super-enhancer để kích hoạt các oncogene — và việc phá vỡ super-enhancer bằng các chất ức chế bromodomain như BET inhibitor đang được nghiên cứu lâm sàng tích cực.

Multi-omics integration cho phép xây dựng **tumor-specific GRN**: từ bulk RNA-seq + ChIP-seq + ATAC-seq + Hi-C của mẫu ung thư, ta có thể xác định TF nào đang "lái" khối u đó, từ đó tìm thuốc nhắm vào TF hoặc CRE của nó.

### 7.2. Drug Target Discovery

**Drug target discovery** (khám phá đích thuốc: quá trình xác định phân tử sinh học — thường là protein — mà nếu điều chỉnh hoạt động của nó sẽ có lợi ích điều trị cho một bệnh cụ thể) từ GRN đi theo logic: nếu một TF điều hòa biểu hiện của một lớn các gene bệnh, TF đó là ứng viên drug target có tính đòn bẩy cao.

**Network medicine** (y học mạng lưới: ứng dụng lý thuyết mạng lưới vào sinh học và y học, xem bệnh là sự rối loạn của các module trong mạng lưới sinh học thay vì rối loạn một phân tử đơn lẻ) tận dụng cấu trúc của GRN để **xếp hạng ứng viên drug target** theo độ trung tâm mạng (network centrality), mức độ đặc hiệu với tế bào bệnh, và khả năng druggable của protein. Các phân tích này giúp rút ngắn giai đoạn tìm kiếm đích từ vài năm xuống vài tháng.

### 7.3. Cell Reprogramming Và Tế Bào Gốc

Khám phá của Yamanaka (2006) rằng chỉ cần bốn TF (OCT4, SOX2, KLF4, MYC) để đưa tế bào trưởng thành về trạng thái tế bào gốc đa năng cảm ứng (**iPSC**) là minh chứng hùng hồn nhất về sức mạnh của TF trong việc định hình GRN. Hiểu đầy đủ cơ chế GRN của quá trình reprogramming — TF nào làm nhiệm vụ gì, theo thứ tự nào, qua những CRE nào — mở ra khả năng hướng dẫn tế bào biệt hóa trực tiếp thành loại tế bào mong muốn mà không cần qua giai đoạn tế bào gốc.

Năm 2024–2026, CellOracle, Geneformer, và các framework in silico perturbation đang được dùng để **thiết kế protocol reprogramming mới** tối ưu hơn, giảm số lượng TF cần thiết và tăng hiệu quả chuyển đổi tế bào. Đây là nền tảng cho liệu pháp tế bào thế hệ tiếp theo.

---

## 8. Thách Thức Và Xu Hướng Năm 2026

### 8.1. Các Thách Thức Còn Tồn Tại

**Độ phủ không đồng đều**: GRN được nghiên cứu kỹ nhất ở tế bào người và chuột, chủ yếu trong bối cảnh ung thư và phát triển phôi. GRN của nhiều mô đặc biệt (tim, thận, hệ thần kinh ngoại vi) vẫn còn ít dữ liệu. Atlas single-cell toàn diện đang lấp dần khoảng trống này.

**Dynamics**: hầu hết dữ liệu omics là **snapshot** — ảnh chụp tại một thời điểm. GRN thực sự là một hệ thống động lực học thay đổi theo thời gian. **Metabolic labeling** kết hợp với scRNA-seq (SLAM-seq, scEU-seq) đang mở ra khả năng đo **RNA dynamics** — tốc độ tổng hợp và phân hủy RNA — để suy luận GRN theo chiều thời gian.

**Causality vs. correlation**: phần lớn GRN hiện tại là **associative** — dựa trên tương quan. **Perturbation experiments** (Perturb-seq: kết hợp CRISPR perturbation với scRNA-seq để đo hệ quả biểu hiện gene của việc knock out từng gene) là hướng đi để thiết lập quan hệ **nhân quả** thay vì chỉ tương quan. Thư viện Perturb-seq quy mô genome đang được xây dựng ở nhiều loại tế bào.

### 8.2. Xu Hướng Năm 2026

- **Perturbation atlas**: tổng hợp kết quả Perturb-seq và CRISPR screen quy mô lớn thành cơ sở dữ liệu GRN nhân quả cho nhiều loại tế bào. **REmap** và **ENCODE4** là hai dự án đang dẫn đầu.
- **Multi-modal foundation models**: kết hợp sequence-to-function model (như Enformer/Borzoi) với single-cell foundation model (như scGPT) để mô hình hóa GRN từ DNA sequence đến phenotype tế bào trong một end-to-end framework.
- **Spatial GRN**: tích hợp spatial transcriptomics + spatial ATAC-seq để xây dựng GRN có độ phân giải không gian, áp dụng trực tiếp vào phân tích mô bệnh học.
- **Epigenetic clocks và GRN aging**: liên kết thay đổi GRN theo tuổi tác với **epigenetic aging clock** (đồng hồ biểu sinh học lão hóa: mô hình dự đoán tuổi sinh học dựa trên mẫu methylation DNA) để hiểu và can thiệp quá trình lão hóa ở cấp độ mạng lưới.

---

## Kết Luận

Gene regulatory network không phải là một đặc điểm thêm thắt của sinh học tế bào — nó **là** sinh học tế bào. Toàn bộ bản sắc của một tế bào, khả năng phản ứng với môi trường, và cách bệnh phát sinh đều nằm trong kiến trúc của mạng lưới này. Từ mô hình TF–DNA binding hai chiều đơn giản, GRN nghiên cứu hiện đại đã tiến đến hệ thống multi-omics tích hợp — đo đồng thời chromatin, biểu hiện gene, cấu trúc 3D, và trình tự DNA — và học từ dữ liệu đó bằng foundation model quy mô lớn.

Năm 2026, ranh giới của lĩnh vực này đang dịch chuyển từ câu hỏi "GRN trông như thế nào?" sang câu hỏi thực dụng hơn: "Chúng ta có thể **thiết kế lại** GRN một cách có kiểm soát để điều trị bệnh và tái lập trình tế bào không?" Đó là câu hỏi sẽ định hình y học phân tử thập kỷ tới.

---

## Tài Liệu Tham Khảo

Lambert, S. A., Jolma, A., Campitelli, L. F., Das, P. K., Yin, Y., Albu, M., Chen, X., Taipale, J., Hughes, T. R., & Weirauch, M. T. (2018). The human transcription factors.
*Cell*, *172*(4), 650–665. https://doi.org/10.1016/j.cell.2018.01.029

Buenrostro, J. D., Giresi, P. G., Zaba, L. C., Chang, H. Y., & Greenleaf, W. J. (2013). Transposition of native chromatin for multimodal regulatory profiling and sequencing.
*Nature Methods*, *10*(12), 1213–1218. https://doi.org/10.1038/nmeth.2688

Lieberman-Aiden, E., van Berkum, N. L., Williams, L., Imakaev, M., Ragoczy, T., Telling, A., Amit, I., Lajoie, B. R., Sabo, P. J., Dorschner, M. O., Sandstrom, R., Bernstein, B., Bender, M. A., Groudine, M., Gnirke, A., Stamatoyannopoulos, J., Mirny, L. A., Lander, E. S., & Dekker, J. (2009). Comprehensive mapping of long-range interactions reveals folding principles of the human genome.
*Science*, *326*(5950), 289–293. https://doi.org/10.1126/science.1181369

Aibar, S., González-Blas, C. B., Moerman, T., Huynh-Thu, V. A., Imrichova, H., Hulselmans, G., Rambow, F., Marine, J. C., Geurts, P., Aerts, J., van den Oord, J., Atak, Z. K., Wouters, J., & Aerts, S. (2017). SCENIC: Single-cell regulatory network inference and clustering.
*Nature Methods*, *14*(11), 1083–1086. https://doi.org/10.1038/nmeth.4463

Kamimoto, K., Stringa, B., Hoffmann, C. M., Jindal, K., Solnica-Krezel, L., & Morris, S. A. (2023). Dissecting cell identity via network inference and in silico gene perturbation.
*Nature*, *614*(7949), 742–751. https://doi.org/10.1038/s41586-022-05688-9

Avsec, Ž., Agarwal, V., Visentin, D., Ledsam, J. R., Grabska-Barwinska, A., Taylor, K. R., Assael, Y., Jumper, J., Kohli, P., & Kelley, D. R. (2021). Effective gene expression prediction from sequence by integrating long-range interactions.
*Nature Methods*, *18*(10), 1196–1203. https://doi.org/10.1038/s41592-021-01252-x

Theodoris, C. V., Xiao, L., Chopra, A., Chaffin, M. D., Al Sayed, Z. R., Hill, M. C., Mantineo, H., Brydon, E. M., Zeng, Z., Liu, X. S., & Ellinor, P. T. (2023). Transfer learning enables predictions in network biology.
*Nature*, *618*(7965), 616–624. https://doi.org/10.1038/s41586-023-06139-9

Mimitou, E. P., Cheng, A., Montalbano, A., Haro-Mora, J. J., Legut, M., Robson, P., Bhanu, N. V., Garcia, B. A., Sanjana, N. E., & Satija, R. (2021). Scalable, multimodal profiling of chromatin accessibility, gene expression and protein levels in single cells.
*Nature Biotechnology*, *39*(10), 1246–1258. https://doi.org/10.1038/s41587-021-00927-2
