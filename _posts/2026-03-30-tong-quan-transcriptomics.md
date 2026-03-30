---
layout: post
title:  "Tổng Quan Về Transcriptomics: Lịch Sử, Công Nghệ Và Sự Tiến Hóa"
author: dat
categories: [ Bioinformatics, Genomics ]
image: assets/my_figs/ds/transcriptomics-overview.png
tags: [ theory, genomics, concept, overview, multi-omics ]
---

Hệ phiên mã (**transcriptome**) là tập hợp toàn bộ các phân tử RNA được tạo ra từ hệ gen của một tế bào hoặc một quần thể tế bào tại một thời điểm nhất định. Không giống như DNA gần như cố định trong mọi tế bào của cơ thể, transcriptome liên tục thay đổi để phản ứng với môi trường, chu kỳ tế bào, độ tuổi và tình trạng bệnh lý. Ngành học nghiên cứu hệ phiên mã, được gọi là **Transcriptomics**, đã trở thành một trong những lĩnh vực mũi nhọn nhất của sinh học phân tử hiện đại. Bài viết này sẽ cung cấp một cái nhìn tổng quan về các công nghệ cốt lõi trong transcriptomics, dòng lịch sử phát triển, những câu hỏi sinh học mà chúng giải quyết và cách các công nghệ này đã biến đổi theo thời gian.

## 1. Dòng Thời Gian Phát Triển Của Transcriptomics

Sự phát triển của transcriptomics được đánh dấu bởi những bước nhảy vọt về công nghệ, đi từ việc đo lường một hoặc vài gen đến hàng chục nghìn gen cùng lúc, và từ mức độ quần thể mô cho đến từng tế bào đơn lẻ (xem Hình 1.1).

### 1.1. Thời Kỳ Khởi Thuỷ: Từ Northern Blot Đến RT-qPCR
Vào cuối những năm 1970, **Northern blot** ra đời, cho phép các nhà khoa học phát hiện và định lượng các đoạn RNA chuyên biệt bằng cách sử dụng các đầu dò (probe) đánh dấu. Tuy nhiên, phương pháp này tốn nhiều thời gian và lượng mẫu lớn. Đến dải những năm 1990, **RT-qPCR** (Reverse Transcription Quantitative PCR) xuất hiện, mang lại độ nhạy và độ chính xác cao hơn rất nhiều, nhưng vẫn bị giới hạn ở khả năng phân tích một vài gen cùng lúc. Những công nghệ này giải quyết câu hỏi: *"Liệu một gen cụ thể có được biểu hiện hay không, và biểu hiện nhiều hay ít?"*

### 1.2. Kỷ Nguyên Đa Nhiệm: Microarray
Năm 1995 đánh dấu một bước ngoặt khi công nghệ **Microarray** (vi mảng) được giới thiệu. Microarray sử dụng hàng nghìn đầu dò DNA cố định trên một con chip nhỏ bé. Khi RNA từ mẫu (sau khi chuyển thành cDNA và đánh dấu huỳnh quang) lai với chip, tín hiệu phát sáng sẽ cho biết mức độ biểu hiện. Microarray đã trả lời được câu hỏi: *"Trong hàng nghìn gen được biết trước, gen nào đang tăng hoặc giảm biểu hiện dưới các điều kiện khác nhau?"*. Dù vậy, Microarray bị giới hạn bởi việc chỉ có thể phát hiện các trình tự đã được thiết kế sẵn (known sequences) và có khoảng động học (dynamic range) hẹp do giới hạn về phổ huỳnh quang.

![Hình 1.1. Dòng thời gian phát triển các công nghệ Transcriptomics.](/assets/my_figs/ds/transcriptomics_timeline.png)
*Hình 1.1. Dòng thời gian tiến hóa của các kỹ thuật phân tích hệ phiên mã.*

## 2. Kỷ Nguyên RNA-Seq Và Sự Tiến Hóa Của Các Nền Tảng

Sự xuất hiện của công nghệ giải trình tự thế hệ mới (**NGS** - Next Generation Sequencing) vào giữa những năm 2000 đã thực sự tạo ra một cuộc cách mạng.

### 2.1. Bulk RNA Sequencing (Giải Trình Tự RNA Khối)
**Bulk RNA-seq** (thường gọi đơn giản là RNA-seq) ra mắt vào năm 2008. Khác với Microarray, RNA-seq đọc trực tiếp từng đoạn RNA, cho phép không chỉ đo lường chính xác mức độ biểu hiện gen mà còn phát hiện các bản phiên mã mới, các biến thể cắt nối (alternative splicing), và đột biến (xem Hình 2.1). RNA-seq giải quyết câu hỏi: *"Toàn bộ hệ phiên mã (kể cả những vùng chưa biết) đang hoạt động như thế nào?"*. Bulk RNA-seq lấy trung bình tín hiệu biểu hiện từ hàng ngàn đến hàng triệu tế bào, giúp ta có bức tranh tổng thể vô cùng rõ nét về một mô hoặc khối u.

![Hình 2.1. So sánh nguyên lý độ phân giải giữa Bulk RNA-seq và Single-cell RNA-seq.](/assets/my_figs/ds/rnaseq_vs_scrnaseq.png)
*Hình 2.1. Minh họa sự khác biệt giữa đo lường tín hiệu trung bình (Bulk RNA-seq) và đo lường từng tế bào riêng biệt (scRNA-seq).*

### 2.2. Single-Cell RNA Sequencing (scRNA-seq)
Năm 2009, bài báo đầu tiên về **Single-cell RNA-seq** (Giải trình tự RNA tế bào đơn) được công bố, mở ra một chiều không gian phân tích mới. Khi Bulk RNA-seq đo lường sinh tố như một "cốc sinh tố" đã xay nhuyễn, scRNA-seq cho phép ta nhìn rõ xem trong đó có bao nhiêu quả dâu tây, bao nhiêu quả chuối. 

scRNA-seq sử dụng hạt vi lỏng (microfluidics) hoặc giếng nano để cô lập từng tế bào, gắn barcode duy nhất cho từng tế bào và từng phân tử RNA. Nó giải quyết triệt để câu hỏi: *"Quần thể tế bào này có đồng nhất không? Có bao nhiêu loại tế bào khác nhau tồn tại trong mẫu, và trạng thái biểu hiện của chúng khác nhau ra sao?"*. Điều này đặc biệt ý nghĩa trong nghiên cứu ung thư và miễn dịch học, nơi các tế bào hiếm đóng vai trò quyết định.

## 3. Chân Trời Mới: Spatial Transcriptomics

Nếu scRNA-seq cho chúng ta biết "CÓ NHỮNG AI" trong mô, thì nó lại làm mất đi thông tin về vị trí "HỌ Ở ĐÂU" do quá trình phân tách tế bào làm vỡ cấu trúc mô.

### 3.1. Sự Tích Hợp Không Gian
**Spatial Transcriptomics** (Hệ phiên mã không gian), được xướng tên là Phương pháp của năm vào 2020 bởi tạp chí *Nature Methods*, là thế hệ tiếp theo. Bằng cách sử dụng các slide mang đầu dò gắn barcode không gian (spatial barcode) trực tiếp lên lát cắt mô, chúng ta có thể định vị chính xác vị trí biểu hiện của từng gen. Nó giúp trả lời câu hỏi tĩnh vi mô: *"Các loại tế bào phân bố ở đâu, tương tác với ai, và vi môi trường ảnh hưởng như thế nào đến biểu hiện gen học?"*.

## 4. Bức Tranh Tổng Thể Và Xu Hướng Tương Lai

Sự tiến hóa từ Northern Blot $\rightarrow$ Microarray $\rightarrow$ Bulk RNA-seq $\rightarrow$ scRNA-seq $\rightarrow$ Spatial Transcriptomics phản ánh sự thay đổi về bản chất câu hỏi sinh học: từ việc theo dõi vài gen đơn lẻ đến phân tích quy mô hệ thống, và từ mức độ trung bình quần thể đến việc bảo toàn chi tiết không gian từng tế bào.

Mỗi công nghệ mới không thay thế hoàn toàn công nghệ gần nhất; thay vào đó, chúng bổ sung cho nhau. Bulk RNA-seq vẫn là tiêu chuẩn vàng về độ sâu giải trình tự và chi phí cho các phân tích bệnh lý mẫu lớn. Trong khi đó, scRNA-seq và Spatial Transcriptomics đang giải quyết những ẩn số vô hạn về sự phát triển mô và vi môi trường khối u. Tương lai của transcriptomics có thể nhắm tới việc đẩy mạnh phân tích **Long-read RNA-seq** để xác định hoàn chỉnh các đồng phân phiên mã (isoforms) một cách chính xác nhất.

## Kết Luận
Transcriptomics đã có một hành trình dài và ấn tượng. Qua nhiều thập kỷ, chúng ta luôn đi tìm các công nghệ có độ phân giải cao hơn, khối lượng dữ liệu khổng lồ hơn để lật mở bức màn giải phẫu ở cấp độ phân tử. Sự ra đời của Microarray, rồi RNA-seq đến phân tích tế bào đơn và hệ phiên mã không gian đã định hình lại nền y sinh và khoa học dữ liệu hiện đại. Chúng tôi sẽ đề cập quy trình thực hành xử lý dữ liệu RNA-seq chi tiết trong bài tutorial riêng.

## Tài Liệu Tham Khảo

Stark, R., Grzelak, M., & Hadfield, J. (2019). RNA sequencing: The teenage years. *Nature Reviews Genetics*, *20*(11), 631–656. https://doi.org/10.1038/s41576-019-0150-2

Wang, Z., Gerstein, M., & Snyder, M. (2009). RNA-Seq: a revolutionary tool for transcriptomics. *Nature Reviews Genetics*, *10*(1), 57–63. https://doi.org/10.1038/nrg2484

Hwang, B., Lee, J. H., & Bang, D. (2018). Single-cell RNA sequencing technologies and bioinformatics pipelines. *Experimental & Molecular Medicine*, *50*(8), 1-14. https://doi.org/10.1038/s12276-018-0071-8

Liao, J., Lu, X., Shao, X., Zhu, L., & Fan, X. (2021). Uncovering an organ’s molecular architecture at single-cell resolution by spatial transcriptomics. *Trends in Biotechnology*, *39*(1), 43–58. https://doi.org/10.1016/j.tibtech.2020.05.006