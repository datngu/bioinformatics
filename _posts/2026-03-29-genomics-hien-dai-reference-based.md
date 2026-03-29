---
layout: post
title:  "Tư Duy Lớn Của Genomics Hiện Đại: Reference Genome Và Hệ Sinh Thái Đa Omics"
author: dat
categories: [ Bioinformatics, Genomics ]
image: assets/my_figs/ds/genomics_hien_dai.png
tags: [ genomics, theory, concept, multi-omics, overview, sequencing ]
---

Nếu phải chọn một ý tưởng duy nhất làm nền móng cho toàn bộ sinh học phân tử hiện đại, đó chính là ý tưởng về **bộ gen tham chiếu** — reference genome. Không phải một công nghệ cụ thể, không phải một phần mềm hay một công trình nghiên cứu đơn lẻ: đó là một **hệ tọa độ dùng chung** mà nhờ nó, mọi câu hỏi về ADN, RNA, protein hay cấu trúc nhiễm sắc thể đều có thể được đặt trong cùng một ngôn ngữ. Bài viết này sẽ giúp bạn hiểu tư duy này từ gốc rễ, và nhìn thấy tại sao hàng chục công nghệ khác nhau — tưởng như rời rạc — thực ra là những mảnh ghép của cùng một bức tranh.

## 1. Reference Genome

Trước khi Dự án Bộ Gen Người (Human Genome Project) hoàn thành vào năm 2003, mỗi phòng thí nghiệm nghiên cứu một đoạn ADN nhỏ một cách cô lập. Giống như hàng nghìn người mù mỗi người đang cầm một mảnh bản đồ khác nhau, nhưng không ai biết mảnh của mình ứng với vị trí nào trên tổng thể.

Reference genome thay đổi tất cả điều đó. Đây là một **bản đồ hệ gen chuẩn** — một trình tự ADN hoàn chỉnh, được lắp ráp, đánh số vị trí và công bố công khai cho từng loài. Đối với người, phiên bản hiện tại được gọi là GRCh38 (hay hg38), bao gồm khoảng 3,2 tỷ cặp base trải dài trên 23 cặp nhiễm sắc thể.

Giá trị của reference genome không nằm ở chỗ nó "hoàn hảo" hay "đại diện cho mọi người". Giá trị của nó nằm ở chỗ nó là **điểm quy chiếu dùng chung**. Khi một nhà nghiên cứu ở Hà Nội và một nhà nghiên cứu ở Boston cùng nói đến "vị trí 179,148,114 trên nhiễm sắc thể 17", họ biết chính xác mình đang nói về cùng một chỗ trên bản đồ. Đây là nền tảng của mọi hợp tác và so sánh trong sinh học phân tử hiện đại.

![Hệ tọa độ reference genome và các lớp omics]({{ site.baseurl }}/assets/my_figs/ds/genomics_hien_dai.png)
*Reference genome đóng vai trò là hệ trục tọa độ chung, nơi mọi lớp thông tin omics đều được neo vào.*

---

## 2. Tại Sao Cần Nhiều Công Nghệ?

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

## 3. Hai Thế Hệ Giải Trình Tự

Để "đo lường" bất kỳ lớp thông tin nào ở trên, ta đều phải giải câu đố tương tự: làm thế nào đọc được các phân tử sinh học một cách tin cậy và có hệ thống? Câu trả lời phổ biến nhất là **giải trình tự ADN** — chuyển đổi phân tử thành dữ liệu số.

### 3.1. Short-read Sequencing

Công nghệ Illumina, đại diện tiêu biểu của thế hệ này, đọc các đoạn ADN ngắn (100–300 nucleotide mỗi đoạn) nhưng với số lượng cực kỳ lớn — hàng trăm triệu đoạn trong một lần chạy. Triết lý ở đây là **độ sâu bù cho độ dài**: từng đoạn riêng lẻ ngắn, nhưng khi chồng lấp hàng triệu đoạn lên nhau, ta có thể xác định trình tự chính xác với độ tin cậy rất cao.

Điều này tạo ra một điểm mạnh và một điểm yếu rõ rệt. Điểm mạnh: độ chính xác gần như tuyệt đối (~99.9%), chi phí thấp, thông lượng cao. Điểm yếu: đoạn ngắn không thể vượt qua được những vùng ADN lặp đi lặp lại dài (như các vùng centromere), giống như bạn muốn lắp ghép một bức tranh ghép hình nhưng tất cả các mảnh đều trông giống nhau.

### 3.2. Long-read Sequencing

PacBio và Oxford Nanopore đại diện cho triết lý ngược lại: thay vì đọc nhiều đoạn ngắn, hãy đọc ít đoạn hơn nhưng **mỗi đoạn dài hơn nhiều** — có thể lên đến hàng chục nghìn, thậm chí hàng triệu nucleotide.

Hãy nghĩ đến sự khác biệt giữa đọc tin nhắn và đọc tiểu thuyết: short-read giống như nhận một nghìn mảnh vụn của nhiều trang sách khác nhau; long-read giống như đọc từng chương liên tục. Đối với những vùng bộ gen phức tạp, lặp lại, hoặc khi cần hiểu cấu trúc tổng thể, long-read là không thể thay thế.

Một khả năng nữa của long-read, đặc biệt là Oxford Nanopore, là **đọc trực tiếp các tín hiệu hóa học** trên ADN mà không cần biến đổi trước — điều này cho phép phát hiện methylation (một dạng biến đổi epigenetic) ngay trong quá trình giải trình tự, thay vì phải làm thí nghiệm riêng.

Trong thực hành, hai thế hệ này **bổ sung cho nhau** nhiều hơn là cạnh tranh: long-read cung cấp khung xương cấu trúc, short-read cung cấp độ chính xác chi tiết.

---

## 4. Các Lớp Omics

Đây là điểm cốt lõi của toàn bộ bài viết này. Mỗi lớp omics không phải là một thế giới riêng biệt — chúng đều **chia sẻ cùng một hệ tọa độ**: vị trí trên reference genome. Điều này biến chúng từ những mảnh thông tin rời rạc thành những lớp của cùng một thực thể.

### 4.1. Genomics

Mỗi người trong chúng ta khác với reference genome khoảng 4–5 triệu vị trí (trong tổng số 3,2 tỷ). Những khác biệt này gọi là **biến thể di truyền**. Một số biến thể lành tính — chỉ là sự đa dạng bình thường giữa người với người. Một số khác, đặc biệt khi xuất hiện ở tế bào ung thư mà không có trong tế bào lành, là tín hiệu quan trọng về cơ chế bệnh.

Việc so sánh bộ gen của tế bào ung thư với tế bào lành từ cùng một người cho phép nhận diện các **đột biến somatic** — những biến đổi ADN chỉ xuất hiện trong khối u, không di truyền, nhưng là động cơ thúc đẩy tế bào phân chia không kiểm soát.

### 4.2. Transcriptomics

ADN là bản thảo, RNA là thông điệp đang được đọc to. Không phải mọi gen đều hoạt động ở mọi lúc và mọi tế bào — bộ gen có tính **ngữ cảnh**. Tế bào cơ tim biểu hiện một tập hợp gen khác tế bào thần kinh, dù cả hai có cùng bộ ADN.

RNA-seq đo lường **lượng RNA** được tạo ra từ mỗi gen tại một thời điểm cụ thể — đây là ảnh chụp nhanh về trạng thái "biểu hiện" của toàn bộ bộ gen. Khi so sánh tế bào bệnh và tế bào lành, ta có thể thấy gen nào đang bị "bật" quá mức, gen nào đang "tắt" bất thường, và từ đó suy ra con đường sinh học nào đang bị rối loạn.

**Single-cell RNA-seq (scRNA-seq)** đẩy khái niệm này lên một tầm cao mới: thay vì lấy trung bình tín hiệu của hàng triệu tế bào trong một mẫu mô, ta đọc biểu hiện gen của **từng tế bào riêng lẻ**. Điều này lần đầu tiên cho phép chúng ta nhìn thấy sự đa dạng ẩn bên trong một mô — phát hiện ra những tiểu quần thể tế bào hiếm, theo dõi quá trình biệt hóa, hoặc xác định tế bào nào trong khối u đang kháng thuốc.

### 4.3. Epigenomics

ADN không nằm trần trụi trong tế bào. Nó cuộn chặt quanh các protein gọi là **histone**, tạo thành cấu trúc bao gói gọi là chromatin. Trạng thái của lớp bao gói này — mở hay đóng, được đánh dấu hóa học như thế nào — quyết định gen nào có thể được đọc và gen nào bị "khóa" lại.

Đây chính là tầng **epigenomics**: những biến đổi không thay đổi trình tự ADN nhưng thay đổi cách bộ gen được sử dụng. Hai tế bào với cùng một bộ ADN có thể hoạt động hoàn toàn khác nhau vì epigenome của chúng khác nhau.

Các công nghệ như **ATAC-seq** xác định vùng nào của chromatin đang mở (tức là đang hoạt động điều hòa), **ChIP-seq** tìm nơi các protein điều hòa bắt vào ADN, còn **WGBS** đo mức độ methyl hóa (một loại "nhãn dán" hóa học) trên từng cytosine trong toàn bộ bộ gen. Tất cả những tín hiệu này đều được ánh xạ về cùng một tọa độ trên reference genome.

Đặc biệt quan trọng: **epigenome có thể thay đổi theo môi trường và bệnh lý**, trong khi genome thì không (ngoại trừ đột biến). Điều này khiến epigenomics trở thành cầu nối giữa di truyền học và môi trường sống.

### 4.4. Hi-C và Cấu Trúc Không Gian Ba Chiều

Bộ gen không phải là một sợi chỉ thẳng. Trong không gian tế bào, ADN gấp lại thành những cấu trúc ba chiều phức tạp, và vị trí không gian này có vai trò quan trọng trong điều hòa gen. Hai vùng ADN cách xa nhau hàng triệu base pair trên trình tự có thể đứng sát nhau trong không gian và tương tác với nhau để bật/tắt gen.

Công nghệ **Hi-C** chụp ảnh những tương tác này — một bản đồ về ai đang "gặp gỡ" ai trong không gian tế bào nhân. Kết quả là chúng ta hiểu tại sao một đột biến ở vùng này lại ảnh hưởng đến gen ở vùng khác tưởng như "xa" trên trình tự thẳng.

---

## 5. Tích Hợp Đa Omics

Sức mạnh thực sự của genomics hiện đại không đến từ từng lớp omics riêng lẻ, mà từ khả năng **tích hợp nhiều lớp** để trả lời những câu hỏi mà không một lớp nào có thể trả lời một mình.

Hãy lấy ví dụ về ung thư vú: một khối u kháng thuốc điều trị. Với chỉ genomics, ta thấy có đột biến trong gen PIK3CA — nhưng đột biến này tồn tại trong cả tế bào nhạy thuốc và kháng thuốc. Transcriptomics cho thấy con đường PI3K/AKT bị kích hoạt mạnh hơn trong nhóm kháng thuốc — nhưng điều gì đang kích hoạt nó? ATAC-seq tiết lộ một vùng enhancer gần gen kháng thuốc đang "mở" ra trong tế bào kháng thuốc — nhưng tại sao? ChIP-seq xác nhận một yếu tố phiên mã đang bám vào vùng đó. Và scRNA-seq phát hiện rằng chỉ khoảng 5% tế bào trong khối u — một tiểu quần thể ẩn — mang toàn bộ đặc điểm kháng thuốc này.

Không có lớp thông tin nào trong số trên, nếu đứng một mình, có thể vẽ ra bức tranh đầy đủ. Chính nhờ chúng được neo vào cùng một hệ tọa độ reference genome, ta có thể chồng lấp và đối chiếu chúng để tìm ra cơ chế.

---

## 6. Tư Duy Từ Câu Hỏi Sinh Học

Một điều quan trọng cần tránh khi mới bắt đầu với tin sinh học là **tư duy từ công nghệ**: "Tôi có dữ liệu RNA-seq, tôi nên làm gì với nó?" Tư duy đúng đắn hơn là **từ câu hỏi sinh học**: "Tôi muốn hiểu tại sao những tế bào này hoạt động khác nhau — dữ liệu nào sẽ giúp tôi trả lời điều đó?"

Câu hỏi sinh học quyết định lớp thông tin cần thiết. Lớp thông tin quyết định công nghệ. Công nghệ quyết định quy trình phân tích. Đây là chuỗi tư duy nên đi từ trên xuống, không phải từ dưới lên.

Reference genome là nền tảng bất biến trong chuỗi này. Dù bạn đang nghiên cứu bộ gen của người, chuột, lúa hay vi khuẩn đường ruột — nguyên lý đều như nhau: xây dựng bản đồ tham chiếu, rồi đo lường mọi thứ **tương đối so với bản đồ đó**.

---

## Kết Luận

Genomics hiện đại không phải là một tập hợp công nghệ ngẫu nhiên. Đó là một **kiến trúc tư duy** được xây dựng xoay quanh một trung tâm duy nhất: reference genome. Mỗi công nghệ omics — từ giải trình tự DNA đến đo trạng thái chromatin, từ đọc RNA đến vẽ bản đồ cấu trúc 3D — là một cách nhìn vào cùng một thực thể từ một góc độ khác nhau. Và nhờ tất cả đều dùng chung một hệ tọa độ, ta có thể tích hợp chúng lại để hiểu những câu hỏi phức tạp nhất của sinh học — từ cơ chế ung thư đến tiến hóa, từ phát triển phôi đến bệnh hiếm gặp.

Trong các bài viết tiếp theo, chúng ta sẽ đi sâu vào từng lớp omics riêng biệt: cơ chế hoạt động của RNA-seq, logic đằng sau ATAC-seq, và triết lý của single-cell genomics. Những bài viết hướng dẫn thực hành (tutorial) sẽ đồng hành để giúp bạn tự tay phân tích dữ liệu thực tế.

## Tài Liệu Tham Khảo

Claussnitzer, M., Cho, J. H., Collins, R., Cox, N. J., Dermitzakis, E. T., Hurles, M. E., … McCarthy, M. I. (2020). A brief history of human disease genetics. *Nature*, *577*(7789), 179–189. https://doi.org/10.1038/s41586-019-1879-7

Stark, R., Grzelak, M., & Hadfield, J. (2019). RNA sequencing: The teenage years. *Nature Reviews Genetics*, *20*(11), 631–656. https://doi.org/10.1038/s41576-019-0150-2

Buenrostro, J. D., Wu, B., Chang, H. Y., & Greenleaf, W. J. (2015). ATAC-seq: A method for assaying chromatin accessibility genome-wide. *Current Protocols in Molecular Biology*, *109*(1), 21.29.1–21.29.9. https://doi.org/10.1002/0471142727.mb2129s109

Lieberman-Aiden, E., van Berkum, N. L., Williams, L., Imakaev, M., Ragoczy, T., Telling, A., … Dekker, J. (2009). Comprehensive mapping of long-range interactions reveals folding principles of the human genome. *Science*, *326*(5950), 289–293. https://doi.org/10.1126/science.1181369

ENCODE Project Consortium. (2012). An integrated encyclopedia of DNA elements in the human genome. *Nature*, *489*(7414), 57–74. https://doi.org/10.1038/nature11247
