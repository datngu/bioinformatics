---
layout: post
title: "Từ Bộ Gen đến Thuốc Đích: Vai Trò Của Omics Trong Nghiên Cứu Ung Thư"
author: dat
categories: [ Bioinformatics, Genomics ]
image: assets/my_figs/ds/omics-ung-thu-drug-target.png
tags: [ genomics, theory, concept, overview, cancer, multi-omics, transcriptomics, epigenomics ]
---

Mỗi năm, ung thư cướp đi gần 10 triệu sinh mạng trên toàn thế giới. Điều đáng chú ý không phải chỉ là con số đó, mà là chúng ta vẫn chưa hiểu đầy đủ tại sao một loại thuốc hiệu quả với bệnh nhân này lại thất bại hoàn toàn với một bệnh nhân khác mang cùng "tên bệnh". Câu trả lời nằm ở một nhận thức mang tính cách mạng: **ung thư không phải là một bệnh duy nhất**. Ung thư là hàng trăm bệnh phân tử khác nhau, và để điều trị hiệu quả từng loại, ta cần công cụ đủ tinh vi để nhìn vào bên trong tế bào ở nhiều chiều thông tin cùng lúc. Đó là sứ mệnh của **nghiên cứu omics** (omics research: nghiên cứu hệ thống và toàn diện một lớp phân tử sinh học cụ thể trong tế bào hoặc mô).

## 1. Bài Toán Ung Thư

### 1.1. Tại Sao Ung Thư Khó Điều Trị

Ung thư phát sinh khi tế bào tích lũy đủ nhiều **đột biến soma** (somatic mutation: biến đổi DNA xuất hiện trong tế bào thân trong suốt cuộc đời, không di truyền sang thế hệ con) để phá vỡ các cơ chế kiểm soát tăng trưởng bình thường. Một khối u điển hình mang hàng nghìn đột biến soma, nhưng phần lớn trong số đó là **đột biến hành khách** (passenger mutation: đột biến xảy ra ngẫu nhiên trong quá trình phân bào nhưng không đóng góp trực tiếp vào sự phát triển ung thư). Chỉ một số ít là **đột biến thúc đẩy** (driver mutation: đột biến tạo ra lợi thế sinh trưởng và chính là động lực thực sự của quá trình sinh ung thư). Nhiệm vụ trung tâm của genomics ung thư là tách bạch hai loại đột biến này.

Khó khăn thứ hai là **tính dị chất nội khối u** (intratumor heterogeneity: ITH, hiện tượng các tế bào trong cùng một khối u không mang cùng bộ đột biến). Một khối u phổi khi sinh thiết có thể chứa nhiều **dòng vô tính** (clones: các quần thể tế bào ung thư đều xuất phát từ một tổ tiên chung nhưng đã tích lũy thêm các đột biến riêng biệt trong quá trình tiến hóa). Khi bệnh nhân điều trị bằng thuốc nhắm trúng đích, dòng vô tính kháng thuốc thiểu số có thể không bị tiêu diệt và dần phát triển thành khối u tái phát — đây là lý do tại sao đơn trị liệu thất bại với nhiều loại ung thư tiến triển.

Hanahan và Weinberg (2011) đã hệ thống hóa những biến đổi căn bản mà tế bào ung thư cần đạt được, gọi là **các đặc trưng ung thư** (hallmarks of cancer): tự duy trì tín hiệu tăng trưởng, kháng các tín hiệu ức chế tăng trưởng, chống lại cái chết tế bào, đạt khả năng phân bào vô hạn, kích thích hình thành tân mạch, xâm lấn và di căn, tái lập trình trao đổi chất, và né tránh hệ miễn dịch. Mỗi đặc trưng này có thể được đạt đến bằng nhiều con đường phân tử khác nhau — và các công nghệ omics lần lượt tiết lộ từng con đường đó.

### 1.2. Giới Hạn Của Y Học Truyền Thống

Trước kỷ nguyên omics, chẩn đoán ung thư dựa chủ yếu vào **giải phẫu bệnh học** (histopathology: quan sát hình thái tế bào dưới kính hiển vi sau khi nhuộm mô bằng hematoxylin-eosin). Bác sĩ giải phẫu bệnh quan sát hình dạng nhân tế bào, mức độ biệt hóa mô, và kiến trúc mô học để phân loại khối u. Kỹ thuật này có giá trị không thể phủ nhận nhưng mang hai giới hạn căn bản.

Thứ nhất, nó chỉ thấy được những gì kính hiển vi quang học có thể thấy — hình dạng và cấu trúc — nhưng không thể cho biết gen nào đang kích hoạt con đường tín hiệu nào, hay protein nào đang giúp tế bào ung thư né tránh hệ miễn dịch.

Thứ hai, và quan trọng hơn, nó **giả định trước** rằng những thứ trông giống nhau thì hành xử giống nhau. Điều này sai. Hai khối u vú có hình thái học giống hệt nhau dưới kính hiển vi có thể có cơ chế phân tử hoàn toàn khác nhau, với tiên lượng và đáp ứng điều trị trái chiều hoàn toàn. **Hóa mô miễn dịch** (immunohistochemistry, IHC: dùng kháng thể đánh dấu để phát hiện protein cụ thể trong mô) cải thiện tình hình bằng cách mở rộng sang một số dấu ấn phân tử, nhưng vẫn bị giới hạn bởi số lượng kháng thể có sẵn và đòi hỏi nhà nghiên cứu phải biết trước muốn tìm gì.

Omics thay đổi quy tắc trò chơi bằng cách cho phép **đo lường toàn hệ thống mà không cần giả thuyết trước** — đây chính là điều mà không kỹ thuật sinh học phân tử cổ điển nào có thể làm được.

---

## 2. Omics và Bộ Ngôn Ngữ Phân Tử Ung Thư

### 2.1. Các Chiều Thông Tin Sinh Học

Nguyên tắc trung tâm của sinh học phân tử mô tả dòng chảy thông tin theo chiều: **DNA → RNA → Protein → Chức năng tế bào**. Mỗi bước chuyển hóa này là một lớp thông tin riêng biệt, và trong ung thư, mỗi lớp có thể bị phá vỡ theo cách riêng của nó. **Omics** là tập hợp các ngành khoa học đo lường toàn hệ thống từng lớp thông tin đó — không đo một hoặc vài phân tử được chọn trước mà đo **toàn bộ** DNA, RNA, protein, hay chất chuyển hóa của một mẫu sinh học. Kết hợp nhiều lớp omics lại gọi là **đa omics** (multi-omics).

| Câu hỏi cần trả lời | Lớp đo lường | Ngành |
|---------------------|-------------|-------|
| DNA ung thư thay đổi những gì? | Biến thể DNA | Genomics |
| Gen nào đang được bật và tắt? | Lượng RNA | Transcriptomics |
| Vùng DNA nào bị "khóa im lặng"? | Trạng thái epigenome | Epigenomics |
| Protein nào đang thực sự hoạt động? | Lượng và trạng thái protein | Proteomics |
| Tế bào ung thư dùng năng lượng thế nào? | Chất chuyển hóa | Metabolomics |

*Bảng 2.1. Ánh xạ từ câu hỏi sinh học ung thư đến lớp omics tương ứng.*

### 2.2. Kiến Trúc Dữ Liệu Chung

Điều làm cho đa omics có thể tích hợp về mặt kỹ thuật là tất cả các lớp dữ liệu đều được định vị trên cùng một **hệ tọa độ reference genome** — mỗi đột biến DNA, mỗi vùng methyl hóa, mỗi tín hiệu biểu hiện gen, mỗi peak chromatin mở đều được neo vào tọa độ nhiễm sắc thể cụ thể. Sự hội tụ tín hiệu từ nhiều lớp omics tại cùng một gen hay con đường tín hiệu mang ý nghĩa xác nhận cao hơn bất kỳ lớp đơn lẻ nào (xem Hình 1).

![Quy trình đa omics trong nghiên cứu ung thư: từ mẫu sinh thiết đến đích thuốc]({{ site.baseurl }}/assets/my_figs/ds/omics-ung-thu-drug-target.png)

*Hình 1. Quy trình nghiên cứu đa omics trong ung thư: từ mẫu sinh thiết qua các lớp omics đến phát hiện đích thuốc.*

---

## 3. Genomics

### 3.1. Phổ Đột Biến Soma

Khi giải trình tự genome của một khối u bằng **whole-genome sequencing** (WGS: giải trình tự toàn genomeome) hay **whole-exome sequencing** (WES: giải trình tự chỉ các vùng mã hóa protein, tức exome, chiếm khoảng 1,5% genome nhưng chứa phần lớn đột biến chức năng), ta thu được một bức tranh toàn diện về các biến đổi DNA soma. Phổ biến nhất là **biến thể nucleotide đơn** (single nucleotide variant, SNV: thay đổi một nucleotide), **indel** (insertion/deletion: thêm hoặc bớt một đến vài nucleotide), và **biến thể số bản sao** (copy number variation, CNV: một đoạn DNA lớn bị nhân lên nhiều lần hay bị xóa đi).

Đặc biệt quan trọng là **biến thể cấu trúc** (structural variant, SV: tái sắp xếp lớn của nhiễm sắc thể như đảo đoạn, chuyển đoạn, hay sáp nhập đoạn lớn). Loại biến thể này thường tạo ra **gen dung hợp** (fusion gene: hai gen riêng biệt bị ghép lại thành một gen mới, tạo ra protein lai có hoạt tính bất thường). Ví dụ điển hình nhất trong lịch sử ung thư học phân tử là sự chuyển đoạn cân bằng giữa nhiễm sắc thể 9 và 22 — tạo ra **nhiễm sắc thể Philadelphia** — trong bệnh bạch cầu mạn dòng tủy (CML, chronic myeloid leukemia). Chuyển đoạn này dẫn đến gen dung hợp **BCR-ABL** mã hóa cho một tyrosine kinase hoạt động liên tục mà không có cơ chế kiểm soát, đây là đột biến thúc đẩy trung tâm của CML.

### 3.2. Từ Đột Biến Driver đến Đích Thuốc

Xác định driver mutation không chỉ là học thuật — nó trực tiếp chỉ ra **đích thuốc** (drug target: phân tử sinh học mà thuốc cần tác động vào để đạt hiệu quả điều trị). Nếu một protein bất thường do đột biến là yếu tố thiết yếu cho sự sống còn của tế bào ung thư, thì ức chế protein đó sẽ tiêu diệt tế bào ung thư trong khi không ảnh hưởng đáng kể đến tế bào bình thường (vốn không có đột biến đó).

**Imatinib** (Gleevec) — thuốc ức chế BCR-ABL — là bằng chứng đầu tiên và hùng hồn nhất: tỷ lệ đạt đáp ứng huyết học hoàn toàn ở bệnh nhân CML tăng từ dưới 5% với hóa trị truyền thống lên **87%** với Imatinib (Druker et al., 2001). Một bước nhảy vọt chưa từng có trong lịch sử ung thư học — và nó chỉ có thể xảy ra khi genomics xác định được BCR-ABL là driver mutation trung tâm.

Những ví dụ tiếp theo không kém phần ấn tượng: đột biến **EGFR** (Epidermal Growth Factor Receptor) trong ung thư phổi không tế bào nhỏ (NSCLC) dẫn đến thuốc ức chế tyrosine kinase thế hệ ba là Osimertinib; **BRAF V600E** trong u ác tính da (melanoma) dẫn đến Vemurafenib; khuếch đại số bản sao của **HER2** trong ung thư vú dẫn đến kháng thể đơn dòng Trastuzumab (Herceptin); và dung hợp **ALK** trong ung thư phổi dẫn đến Crizotinib. Toàn bộ nhóm thuốc nhắm trúng đích này đều bắt nguồn từ một thực nghiệm genomics đã nhìn thấy bất thường ở cấp độ DNA.

### 3.3. Sinh Thiết Lỏng

Một ứng dụng mà hình thức chẩn đoán truyền thống không thể thực hiện được là **sinh thiết lỏng** (liquid biopsy). Khi tế bào ung thư chết, chúng giải phóng các mảnh DNA vào máu, gọi là **DNA ung thư lưu thông** (circulating tumor DNA, ctDNA). Bằng cách giải trình tự DNA huyết tương với độ nhạy cực cao, ta có thể phát hiện đột biến ung thư từ một ống máu thông thường, theo dõi đáp ứng điều trị theo thời gian thực, và phát hiện sự xuất hiện của đột biến kháng thuốc trước khi ung thư tái phát lâm sàng hàng tháng. Đây là ứng dụng hoàn toàn nằm ngoài tầm với của bất kỳ kỹ thuật chẩn đoán không dựa trên giải trình tự nào.

---

## 4. Transcriptomics

### 4.1. Phân Loại Phân Tử Ung Thư

Trong khi genomics hỏi "DNA thay đổi gì?", transcriptomics hỏi "Genome đang *làm gì*?" Mức biểu hiện gen phản ánh trạng thái hoạt động của tế bào tại thời điểm đo, và điều này mang thông tin lâm sàng mà DNA không thể cung cấp. Nhiều khối u không mang đột biến DNA rõ ràng nhưng có hồ sơ biểu hiện gen rối loạn nghiêm trọng, và điều đó chỉ transcriptomics mới nhìn thấy được.

Minh chứng đột phá nhất đến từ nghiên cứu ung thư vú của Perou và cộng sự (2000): bằng cách phân tích cùng lúc biểu hiện của hàng nghìn gen, họ phát hiện rằng ung thư vú có thể được chia thành **bốn phân nhóm phân tử** — Luminal A, Luminal B, HER2-enriched, và Basal-like — với tiên lượng và đáp ứng điều trị hoàn toàn khác nhau. Hai khối u Luminal A và Basal-like có thể trông giống nhau dưới kính hiển vi nhưng có tỷ lệ sống sót 10 năm chênh lệch đến hơn 30 điểm phần trăm. Điều đó không thể nhận ra bằng bất kỳ kỹ thuật hình thái học nào. Phân loại phân tử này ngày nay được chuẩn hóa thành genome **PAM50** và được sử dụng lâm sàng để quyết định liệu bệnh nhân có cần hóa trị bổ trợ hay không.

### 4.2. Giải Phẫu Khối U Ở Cấp Tế Bào Đơn

**Giải trình tự RNA đơn tế bào** (single-cell RNA sequencing, scRNA-seq: đo biểu hiện gen của từng tế bào riêng lẻ thay vì lấy trung bình của cả quần thể tế bào) là bước tiến tiếp theo. Trong khối u, scRNA-seq tiết lộ điều mà RNA-seq khối (bulk RNA-seq) che khuất: không phải tất cả tế bào ung thư đều thực thi cùng một chương trình gen. Một số tế bào biểu hiện gen liên quan đến phân bào nhanh, một số khác biểu hiện gen liên quan đến di căn, và một quần thể nhỏ — **tế bào gốc ung thư** (cancer stem cells: quần thể có khả năng tự đổi mới và kháng thuốc cao, thường là nguồn gốc của tái phát) — biểu hiện chương trình tế bào gốc giúp chúng sống sót qua hóa trị.

scRNA-seq cũng lập bản đồ **vi môi trường khối u** (tumor microenvironment, TME: tập hợp các tế bào không ung thư bao quanh khối u gồm tế bào miễn dịch, tế bào sợi, và tế bào nội mô mạch máu, có ảnh hưởng lớn đến sự phát triển và kháng thuốc của khối u). Tỷ lệ tế bào T điều hòa so với tế bào T gây độc trong vi môi trường là một trong những yếu tố dự đoán đáp ứng với liệu pháp kiểm soát miễn dịch — và chỉ scRNA-seq mới có thể đo lường đồng thời cả hai quần thể trên cùng một mẫu.

---

## 5. Epigenomics

### 5.1. Khóa Im Lặng Gen Mà Không Xóa DNA

Ung thư không chỉ biến đổi DNA — nó còn thay đổi **cách DNA được đọc** mà không làm thay đổi trình tự nucleotide. Đây là lĩnh vực **epigenomics** (nghiên cứu toàn bộ các biến đổi epigenetic: thay đổi biểu hiện gen không do thay đổi trình tự DNA). Hai cơ chế epigenetic trung tâm là **methyl hóa DNA** (DNA methylation: gắn nhóm methyl lên cytosine tại các vị trí CpG, thường dẫn đến tắt biểu hiện gen) và **biến đổi histone** (histone modification: thay đổi hóa học trên đuôi histone, điều tiết mức độ nén chặt của nhiễm sắc thể và khả năng tiếp cận DNA của bộ máy phiên mã).

Khi promoter của một **gen ức chế khối u** (tumor suppressor gene: gen mã hóa protein ức chế tăng trưởng không kiểm soát) bị hypermethylation (methyl hóa quá mức), gen đó bị "khóa im lặng" mà không cần đột biến DNA thực sự — và tế bào ung thư vẫn đạt được kết quả tương tự như khi gen đó bị xóa. Ví dụ: hypermethylation promoter của **MLH1** (gen sửa chữa lỗi bắt cặp DNA) là cơ chế phổ biến dẫn đến **bất ổn định vi vệ tinh** (microsatellite instability, MSI) trong ung thư đại tràng sporadic — một đặc điểm quan trọng vì MSI-high hiện là chỉ định cho liệu pháp kiểm soát miễn dịch.

### 5.2. Đột Biến Tạo Ra Cảnh Quan Epigenome Mới

Mối liên hệ giữa genomics và epigenomics trở nên đặc biệt sâu sắc ở trường hợp **đột biến IDH1/IDH2** (isocitrate dehydrogenase 1 và 2) trong u não cấp độ thấp và bạch cầu cấp dòng tủy. Đột biến gain-of-function ở hai gen này khiến enzyme tạo ra **2-hydroxyglutarate** (2-HG: một **oncometabolite** — chất chuyển hóa bất thường được tạo ra do đột biến, tích lũy đến nồng độ cao bất thường và gây rối loạn chức năng tế bào). 2-HG ức chế cạnh tranh enzyme **TET** (tham gia demethyl hóa DNA) và **histone KDM** (lysine demethylase), dẫn đến **hypermethylation toàn cầu** ảnh hưởng đến biểu hiện hàng trăm gen — được gọi là kiểu hình G-CIMP (glioma CpG island methylator phenotype). Chỉ từ một đột biến điểm duy nhất trong gen IDH, toàn bộ cảnh quan epigenome của tế bào bị biến đổi — đây là ví dụ rõ ràng tại sao genomics và epigenomics không thể nghiên cứu tách rời nhau.

### 5.3. Epigenome Là Đích Thuốc Đặc Biệt

Vì các biến đổi epigenetic không làm thay đổi trình tự DNA, về lý thuyết chúng **có thể đảo ngược** — đây là lý do tại sao nhóm thuốc epigenetic (DNMT inhibitor, HDAC inhibitor, IDH inhibitor) đang được ứng dụng mạnh mẽ. Ivosidenib (ức chế IDH1) và Enasidenib (ức chế IDH2) đã được FDA chấp thuận cho điều trị AML có đột biến IDH, trực tiếp dựa trên hiểu biết epigenomic từ omics. Điều đặc biệt hơn là **methyl hóa DNA trong huyết tương** đang được khai thác để phát hiện sớm ung thư không triệu chứng: mỗi loại ung thư để lại một "dấu ấn methyl hóa" đặc trưng trong ctDNA, làm nền tảng cho các xét nghiệm tầm soát đa ung thư thế hệ mới.

---

## 6. Proteomics và Metabolomics

### 6.1. Khoảng Cách Từ RNA Đến Chức Năng Thực

Genome và transcriptome nói lên tiềm năng và ý định — nhưng **protein mới là người thực sự thực thi** mọi chức năng tế bào. Có một khoảng cách quan trọng giữa mRNA và protein: không phải tất cả mRNA đều được dịch mã với cùng hiệu suất, và không phải tất cả protein đều có cùng hoạt tính sau khi được tổng hợp.

**Biến đổi sau dịch mã** (post-translational modifications, PTMs: các biến đổi hóa học thêm vào protein sau tổng hợp, như phosphoryl hóa, ubiquitin hóa, hay acetyl hóa) là ngôn ngữ tế bào dùng để bật/tắt hoạt tính protein theo thời gian thực. Phosphoryl hóa là biến đổi quan trọng nhất trong ung thư: khi một **kinase** (enzyme gắn nhóm phosphate vào protein đích, thường kích hoạt tín hiệu tế bào) bị đột biến và hoạt động liên tục—như BCR-ABL hay EGFR đột biến—nó gây phosphoryl hóa bất thường hàng trăm protein xuôi dòng. **Phosphoproteomics** (đo lường trạng thái phosphoryl hóa của toàn bộ proteome) vẽ ra bản đồ "ai đang kích hoạt ai" trong tế bào ung thư — thông tin không thể thu được từ genomics hay transcriptomics.

Dự án **CPTAC** (Clinical Proteomic Tumor Analysis Consortium: dự án tích hợp proteomics lâm sàng với dữ liệu TCGA genomics và transcriptomics) của NIH đã chứng minh rằng mức protein thực sự hiện diện không tương quan cao với mức RNA tương ứng ở nhiều gen — nhấn mạnh sự cần thiết không thể bỏ qua của lớp proteomics trong bức tranh đa omics toàn diện.

### 6.2. Trao Đổi Chất Ung Thư

**Metabolomics** (đo lường toàn bộ **chất chuyển hóa** (metabolites: các phân tử nhỏ là sản phẩm trung gian và cuối cùng của quá trình sinh hóa trong tế bào) tiết lộ rằng tế bào ung thư tái lập trình trao đổi chất theo những cách đặc trưng. Hiện tượng **Warburg effect** (hiệu ứng Warburg: tế bào ung thư ưu tiên con đường glycolysis kỵ khí ngay cả khi có đủ oxy, tạo ra lactate thay vì sử dụng chu trình Krebs hiệu quả hơn) được Otto Warburg phát hiện từ năm 1924 nhưng chỉ đến khi metabolomics hiện đại ra đời, ta mới hiểu toàn bộ chiều sâu của việc tái lập trình trao đổi chất ung thư. Như đã đề cập ở phần epigenomics, đột biến IDH1/IDH2 tạo ra 2-HG, một oncometabolite phát hiện bởi metabolomics và là đích thuốc trực tiếp — ví dụ đẹp nhất về vòng khép kín genomics → metabolomics → drug target.

---

## 7. Tích Hợp Đa Omics và Phát Hiện Đích Thuốc

### 7.1. Logic Của Sự Hội Tụ Tín Hiệu

Không có lớp omics nào đủ sức hoàn chỉnh bức tranh ung thư một mình. Sức mạnh thực sự xuất hiện khi nhiều lớp dữ liệu **hội tụ tín hiệu vào cùng một thực thể sinh học**: khi một gen đồng thời bị đột biến driver (genomics), biểu hiện tăng cao bất thường (transcriptomics), promoter không bị methyl hóa trong khi gen ức chế cùng con đường lại bị im lặng epigenetic (epigenomics), và protein tương ứng bị phosphoryl hóa quá mức (proteomics) — đó là tín hiệu cực mạnh rằng gen này là đích thuốc quan trọng.

### 7.2. TCGA và Bộ Dữ Liệu Đa Omics Quy Mô Lớn

**TCGA** (The Cancer Genome Atlas: Dự án Bản Đồ Bộ Gen Ung Thư) là minh chứng rõ ràng nhất cho sức mạnh của đa omics ở quy mô dân số. Với hơn 11.000 bệnh nhân và hơn 33 loại ung thư, TCGA thu thập đồng thời dữ liệu WES, RNA-seq, methyl hóa genome-wide, CNV, và proteomics cho cùng một mẫu khối u. Phân tích tích hợp toàn bộ tập dữ liệu này đã phát hiện những điểm tương đồng phân tử bất ngờ giữa các loại ung thư về mặt mô học: ung thư bàng quang và một phân nhóm ung thư vú chia sẻ cùng một chương trình phân tử gần giống nhau hơn cả hai phân nhóm ung thư cùng cơ quan (Cancer Genome Atlas Research Network, 2013). Đây là lý do để bắt đầu thiết kế nghiên cứu lâm sàng theo phân loại phân tử thay vì vị trí giải phẫu.

### 7.3. Hai Ví Dụ Điển Hình Về Tích Hợp Đa Omics

Ví dụ đầu tiên là **liệu pháp kiểm soát miễn dịch** (immune checkpoint inhibitors: thuốc phong tỏa phân tử ức chế miễn dịch trên tế bào T như PD-1/PD-L1 hay CTLA-4). Việc dự đoán bệnh nhân nào sẽ đáp ứng đòi hỏi ít nhất ba lớp omics: **gánh nặng đột biến khối u** (tumor mutational burden, TMB: tổng số đột biến soma trên mỗi megabase DNA — khối u có TMB cao tạo ra nhiều **tân kháng nguyên** (neoantigen: peptide bất thường từ protein đột biến, bị hệ miễn dịch nhận diện như "ngoại lai") hơn, thu hút tế bào T nhiều hơn; mức biểu hiện PD-L1 trên tế bào ung thư (xác định bởi transcriptomics hoặc IHC); và thành phần tế bào miễn dịch trong vi môi trường (xác định bởi scRNA-seq). Không có lớp omics đơn lẻ nào đủ để dự đoán đáp ứng — chỉ khi tích hợp cả ba chiều nhìn mới cho độ chính xác lâm sàng có ý nghĩa (Topalian et al., 2012).

Ví dụ thứ hai là **tính gây chết tổng hợp** (synthetic lethality: hiện tượng trong đó sự kết hợp của hai biến đổi mới gây chết tế bào, trong khi mỗi biến đổi đơn lẻ thì không). Đột biến **BRCA1/BRCA2** — xác định bởi genomics — làm tế bào ung thư mất khả năng sửa chữa DNA qua cơ chế tái tổ hợp tương đồng (homologous recombination). Kết hợp điều này với ức chế **PARP** (poly-ADP-ribose polymerase: enzyme tham gia sửa chữa DNA theo cơ chế khác) bằng Olaparib gây tích lũy tổn thương DNA không hồi phục trong tế bào BRCA-deficient, tiêu diệt đặc hiệu ung thư trong khi tế bào bình thường (vẫn còn BRCA lành mạnh) sống sót (Bryant et al., 2005). Toàn bộ logic điều trị này bắt nguồn từ hiểu biết genomics sâu về con đường sửa chữa DNA.

---

## 8. Điều Omics Có Thể Làm Mà Kỹ Thuật Khác Không Thể

### 8.1. Khám Phá Không Cần Giả Thuyết Trước

Mọi kỹ thuật sinh học phân tử cổ điển đều **giả thuyết-hướng dẫn** (hypothesis-driven): Western Blot đo protein bạn đã biết tên, RT-PCR đo gen bạn đã biết trình tự, IHC tìm kháng nguyên bạn đã có kháng thể. Omics đảo ngược hoàn toàn quy trình này: nó đo **toàn bộ** rồi để dữ liệu tự tiết lộ điều quan trọng.

Đây là lý do TCGA phát hiện **IDH1/IDH2** là driver mutation quan trọng trong glioma — không ai biết trước để thiết kế thực nghiệm tìm chúng. Dữ liệu genomics trên hàng nghìn mẫu tự để lộ tín hiệu đột biến lặp lại với tần suất không thể là ngẫu nhiên. Đây là mô hình khám phá **data-driven** mà y học truyền thống không thể làm được.

### 8.2. Phân Loại Lại Bệnh Theo Bản Chất Phân Tử

Omics không chỉ giúp điều trị tốt hơn một bệnh đã được định nghĩa — nó **tái định nghĩa chính bệnh đó**. Ung thư vú không còn là một bệnh mà là ít nhất năm thực thể phân tử. Ung thư phổi NSCLC không còn là một bệnh mà là một tập hợp bệnh được phân loại bởi đột biến driver: EGFR-mutant, ALK-fusion, KRAS-mutant, ROS1-fusion, PD-L1-high... Mỗi nhóm có thuốc đích riêng và phác đồ điều trị riêng. Sự tái phân loại này cứu sống bệnh nhân — và toàn bộ nó đến từ omics.

### 8.3. Giám Sát Tiến Hóa Dòng Vô Tính Theo Thời Gian Thực

Ung thư tiến hóa dưới áp lực chọn lọc của thuốc. Dòng vô tính thiểu số mang đột biến kháng thuốc có thể phát triển thành quần thể chiếm ưu thế sau điều trị. Với ctDNA, ta có thể theo dõi **entropy di truyền khối u** theo thời gian thực qua xét nghiệm máu định kỳ — phát hiện đột biến kháng thuốc xuất hiện sớm hơn vài tháng so với tái phát lâm sàng. Điều này cho phép điều chỉnh phác đồ chủ động thay vì phản ứng sau khi ung thư đã tái phát hoàn toàn — hoàn toàn nằm ngoài tầm với của mọi kỹ thuật không dựa trên giải trình tự DNA.

---

## Kết Luận

Ung thư là một cuộc chiến diễn ra trong thế giới vi mô của phân tử, và để chiến thắng, ta cần nhìn vào thế giới đó ở nhiều chiều cùng một lúc. Genomics tiết lộ bản đồ đột biến và chỉ ra đích thuốc; transcriptomics phản ánh trạng thái hoạt động và cho phép phân loại phân tử; epigenomics hé lộ những gen bị im lặng không phải do đột biến; proteomics xác nhận những gì thực sự xảy ra ở cấp chức năng; metabolomics nắm bắt hậu quả trao đổi chất của toàn bộ những biến đổi trên.

Không lớp nào thay thế được lớp kia — chúng là các mặt phẳng khác nhau của cùng một thực thể phức tạp. Sức mạnh của đa omics không nằm ở khối lượng dữ liệu mà ở sự hiểu biết hệ thống xuất hiện khi chúng được tích hợp. Từ BCR-ABL đến BRCA/PARP, từ PD-L1 đến IDH1/2, mỗi đích thuốc thành công trong y học ung thư hiện đại đều bắt nguồn từ một thực nghiệm omics đã nhìn thấy điều mà kính hiển vi và kháng thể không thể.

Chúng tôi sẽ đề cập đến quy trình phân tích bioinformatics cụ thể cho từng lớp omics — từ variant calling, RNA-seq differential expression, đến phân tích methylation và tích hợp đa omics — trong các bài tutorial riêng.

---

## Tài Liệu Tham Khảo

Druker, B. J., Talpaz, M., Resta, D. J., Peng, B., Buchdunger, E., Ford, J. M., Lydon, N. B., Kantarjian, H., Capdeville, R., Ohno-Jones, S., & Sawyers, C. L. (2001). Efficacy and safety of a specific inhibitor of the BCR-ABL tyrosine kinase in chronic myeloid leukemia.
*New England Journal of Medicine*, *344*(14), 1031–1037. https://doi.org/10.1056/NEJM200104053441401

Hanahan, D., & Weinberg, R. A. (2011). Hallmarks of cancer: The next generation.
*Cell*, *144*(5), 646–674. https://doi.org/10.1016/j.cell.2011.02.013

Cancer Genome Atlas Research Network. (2013). The Cancer Genome Atlas Pan-Cancer analysis project.
*Nature Genetics*, *45*(10), 1113–1120. https://doi.org/10.1038/ng.2764

Perou, C. M., Sørlie, T., Eisen, M. B., van de Rijn, M., Jeffrey, S. S., Rees, C. A., Pollack, J. R., Ross, D. T., Johnsen, H., Akslen, L. A., Fluge, Ø., Pergamenschikov, A., Williams, C., Zhu, S. X., Lønning, P. E., Børresen-Dale, A. L., Brown, P. O., & Botstein, D. (2000). Molecular portraits of human breast tumours.
*Nature*, *406*(6797), 747–752. https://doi.org/10.1038/35021093

Bryant, H. E., Schultz, N., Thomas, H. D., Parker, K. M., Flower, D., Lopez, E., Kyle, S., Meuth, M., Curtin, N. J., & Helleday, T. (2005). Specific killing of BRCA2-deficient tumours with inhibitors of poly(ADP-ribose) polymerase.
*Nature*, *434*(7035), 913–917. https://doi.org/10.1038/nature03443

Topalian, S. L., Hodi, F. S., Brahmer, J. R., Gettinger, S. N., Smith, D. C., McDermott, D. F., Powderly, J. D., Carvajal, R. D., Sosman, J. A., Atkins, M. B., Leming, P. D., Spigel, D. R., Antonia, S. J., Horn, L., Drake, C. G., Pardoll, D. M., Chen, L., Sharfman, W. H., Anders, R. A., & Wolchok, J. D. (2012). Safety, activity, and immune correlates of anti-PD-1 antibody in cancer.
*New England Journal of Medicine*, *366*(26), 2443–2454. https://doi.org/10.1056/NEJMoa1200690

Dang, L., White, D. W., Gross, S., Bennett, B. D., Bittinger, M. A., Driggers, E. M., Fantin, V. R., Jang, H. G., Jin, S., Keenan, M. C., Marks, K. M., Prins, R. M., Ward, P. S., Yen, K. E., Liau, L. M., Rabinowitz, J. D., Cantley, L. C., Thompson, C. B., Vander Heiden, M. G., & Su, S. M. (2009). Cancer-associated IDH1 mutations produce 2-hydroxyglutarate.
*Nature*, *462*(7274), 739–744. https://doi.org/10.1038/nature08617
