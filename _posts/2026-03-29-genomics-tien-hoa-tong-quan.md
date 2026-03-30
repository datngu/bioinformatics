---
layout: post
title:  "Genomics Tiến Hoá: Tổng Quan Lý Thuyết và Phương Pháp Hiện Đại"
author: dat
categories: [ Bioinformatics, Genomics ]
image: assets/my_figs/ds/genomics-tien-hoa-tong-quan.png
tags: [ theory, genomics, overview, intermediate ]
---

Tại sao con người và tinh tinh chia sẻ tới 98,7% genome? Tại sao một số quần thể người có khả năng kháng bệnh sốt rét cao hơn những quần thể khác? Tại sao cá voi lưng gù lại gần gũi hơn về mặt di truyền với hươu cao cổ so với cá mập? Những câu hỏi này, tuy xuất phát từ các góc độ khác nhau, đều có chung một nền tảng giải đáp: **genomics tiến hoá** (evolutionary genomics). Đây là lĩnh vực nghiên cứu cách toàn bộ thông tin di truyền của sinh vật biến đổi, tích lũy và được định hình qua hàng triệu năm tiến hoá. Khác với di truyền học truyền thống chỉ nghiên cứu từng gene riêng lẻ, genomics tiến hoá đặt toàn bộ genome vào bối cảnh thời gian sâu và phả hệ tiến hoá, từ đó mở ra một tầm nhìn toàn diện về sự sống.

## 1. Genomics Tiến Hoá Là Gì

### 1.1. Định Nghĩa và Phạm Vi

**Genomics tiến hoá** (evolutionary genomics) là ngành khoa học nghiên cứu cấu trúc, chức năng và sự biến đổi của genome dưới ánh sáng của lý thuyết tiến hoá. Nói cách khác, đây là sự giao thoa giữa hai lĩnh vực lớn: **sinh học tiến hoá** (evolutionary biology), vốn quan tâm đến nguồn gốc và sự biến đổi của các đặc tính sinh học qua thời gian, và **genomics** (genome học), vốn tập trung vào việc giải mã và phân tích tổng thể genome. Sự kết hợp này cho phép các nhà khoa học đặt câu hỏi không chỉ là "gene nào mã hoá protein X?" mà còn "tại sao gene đó lại tồn tại ở dạng này, và không phải dạng khác, qua hàng trăm triệu năm?"

Phạm vi của genomics tiến hoá rất rộng. Nó bao trùm từ việc tái dựng cây phả hệ tiến hoá của các loài, phân tích sự thay đổi tần số alen trong quần thể, đến việc xác định những vùng genome chịu áp lực chọn lọc mạnh. Lĩnh vực này cũng tiếp thu dữ liệu từ nhiều nguồn khác nhau, bao gồm DNA hiện đại, DNA cổ đại được khai thác từ hóa thạch, và siêu genome (**metagenome**) của cộng đồng vi sinh vật.

### 1.2. Lịch Sử Hình Thành

Nền tảng lý thuyết của genomics tiến hoá bắt đầu từ **thuyết tiến hoá tổng hợp hiện đại** (Modern Synthesis) vào giữa thế kỷ 20, khi di truyền học Mendel được tích hợp với lý thuyết chọn lọc tự nhiên của Darwin. Tuy nhiên, lĩnh vực này chỉ thực sự bùng nổ sau khi công nghệ **giải trình tự thế hệ kế tiếp** (next-generation sequencing, NGS) trở nên phổ biến vào đầu thế kỷ 21. Dự án Hệ Gen Người (Human Genome Project) hoàn thành năm 2003 đã đặt nền móng quan trọng, nhưng chính sự ra đời của giải trình tự thông lượng cao mới cho phép so sánh hàng nghìn genome cùng lúc, biến genomics tiến hoá thành một ngành khoa học dữ liệu thực thụ.

Từ cột mốc năm 2008, khi dự án 1000 Genomes bắt đầu thu thập dữ liệu, việc giải trình tự hàng nghìn cá thể người từ các dân số khác nhau đã cung cấp một nguồn tài nguyên vô giá để nghiên cứu chọn lọc tự nhiên, di cư và cấu trúc quần thể ở cấp độ genome toàn phần.

## 2. Các Khái Niệm Cốt Lõi

### 2.1. Bốn Lực Tiến Hoá

Để hiểu genomics tiến hoá, cần nắm vững bốn lực cơ bản định hình sự biến đổi di truyền theo thời gian.

**Đột biến** (mutation) là nguồn gốc của mọi biến dị di truyền. Mỗi lần DNA được sao chép, các lỗi sao chép nhỏ có thể xảy ra, tạo ra các alen mới trong quần thể. Tốc độ đột biến điểm ở người vào khoảng 1,2 × 10⁻⁸ mỗi cặp bazơ mỗi thế hệ, một con số nhỏ nhưng đủ để tạo ra khoảng 35 đột biến mới trong genome của mỗi cá thể.

**Chọn lọc tự nhiên** (natural selection) là quá trình các alen làm tăng khả năng sinh tồn và sinh sản được truyền lại nhiều hơn qua các thế hệ. Chọn lọc có thể ở nhiều dạng: **chọn lọc định hướng** (directional selection) đẩy quần thể về một đặc điểm cực đoan, **chọn lọc duy trì** (purifying selection) loại bỏ các đột biến có hại, và **chọn lọc cân bằng** (balancing selection) duy trì đa dạng alen trong quần thể. Trong genomics tiến hoá, chọn lọc được phát hiện bằng cách tìm kiếm những vùng genome có tần số alen thay đổi nhanh bất thường hoặc thiếu đa dạng di truyền so với phần còn lại của genome.

**Lệch lạc di truyền** (genetic drift) là sự thay đổi tần số alen do may rủi ngẫu nhiên, không liên quan đến giá trị thích nghi của alen đó. Hiệu ứng này đặc biệt mạnh trong các quần thể nhỏ, và có thể dẫn đến mất đi hoàn toàn các alen dù chúng không có hại gì. **Hiệu ứng cổ chai** (bottleneck effect) là một dạng đặc biệt của lệch lạc di truyền, xảy ra khi một quần thể bị thu hẹp đột ngột rồi mở rộng trở lại, ví dụ điển hình là tổ tiên người hiện đại rời châu Phi khoảng 70.000 năm trước.

**Gene flow** (gene flow: sự di chuyển alen giữa các quần thể thông qua lai giống và di cư) là sự di chuyển của alen giữa các quần thể thông qua lai giống và di cư. Gene flow có thể làm đồng nhất hoá các quần thể về mặt di truyền, hoặc đưa vào những alen thích nghi mới. Trường hợp nổi bật là người hiện đại tiếp nhận các gene liên quan đến kháng lạnh và kháng bệnh từ người Denisovan khi di cư qua châu Á, để lại dấu ấn trong genome các dân tộc Đông Nam Á và châu Đại Dương ngày nay.

### 2.2. Đồng Hồ Phân Tử

**Đồng hồ phân tử** (molecular clock) là nguyên lý cho rằng DNA tích lũy đột biến với tốc độ tương đối ổn định theo thời gian. Từ tốc độ đột biến đã biết, các nhà khoa học có thể ước tính thời điểm phân kỳ tiến hoá giữa các loài hoặc các dòng tộc. Ví dụ, bằng cách so sánh sự khác biệt di truyền giữa người và tinh tinh với tốc độ đột biến đã hiệu chỉnh, các nhà nghiên cứu ước tính tổ tiên chung cuối cùng của hai loài sống cách đây khoảng 5 đến 7 triệu năm.

Cần lưu ý rằng đồng hồ phân tử không hoàn toàn đều đặn; tốc độ đột biến có thể thay đổi theo loài, vùng genome và các yếu tố sinh học khác như tuổi sinh sản và kích thước quần thể. Vì vậy, các phương pháp hiện đại thường sử dụng **đồng hồ phân tử thư giãn** (relaxed molecular clock) cho phép tốc độ tiến hoá khác nhau giữa các nhánh của cây tiến hoá.

### 2.3. Cấu Trúc Quần Thể

**Cấu trúc quần thể** (population structure) mô tả sự phân bố của biến dị di truyền trong và giữa các quần thể. Khi một quần thể bị chia cắt về mặt địa lý hoặc sinh thái, các dòng tộc bắt đầu tích lũy đột biến độc lập với nhau, dẫn đến sự phân hoá di truyền. Mức độ phân hoá này được đo bằng hệ số **F*ST* (fixation index)**, dao động từ 0 (không có sự phân hoá) đến 1 (phân hoá hoàn toàn).

Công cụ phân tích thành phần chính (**Principal Component Analysis**, PCA) và mô hình phân cụm xác suất (điển hình là ADMIXTURE) được dùng rộng rãi để trực quan hoá cấu trúc quần thể, và đã giúp tái dựng lịch sử di cư của loài người một cách chi tiết và chính xác.

## 3. Phương Pháp Hiện Đại

### 3.1. Genomics Quần Thể

**Genomics quần thể** (population genomics) mở rộng di truyền học quần thể truyền thống lên cấp độ toàn genome. Thay vì phân tích một vài locus, các nghiên cứu hiện đại giải trình tự hàng nghìn cá thể và quét hàng triệu **đa hình nucleotide đơn** (single nucleotide polymorphism, SNP) để phát hiện tín hiệu chọn lọc, phân tán quần thể và nguồn gốc phả hệ. Dự án 1000 Genomes, gnomAD và UK Biobank là những ví dụ tiêu biểu cho cách tiếp cận này ở người.

Một khái niệm quan trọng trong genomics quần thể là **cân bằng liên kết** (linkage disequilibrium, LD), chỉ sự liên kết phi ngẫu nhiên giữa các alen ở các locus khác nhau. Khi một alen mới xuất hiện do đột biến, nó mang theo một nền chuỗi haplotype xung quanh. Tốc độ phân rã của LD theo thời gian do tái tổ hợp cho phép ước tính thời điểm chọn lọc xảy ra, và là cơ sở cho các phương pháp phát hiện chọn lọc như iHS và XP-EHH.

### 3.2. Phát Sinh Loài Học Hệ Gen

**Phát sinh loài học genome** (phylogenomics) xây dựng cây tiến hoá dựa trên hàng trăm đến hàng nghìn gene, thay vì một hoặc vài gene như phương pháp cổ điển. Điều này giúp giải quyết những nhánh tiến hoá trước đây không thể phân giải do thiếu tín hiệu di truyền. Nghiên cứu của Jarvis et al. (2014) về phát sinh loài của 48 loài chim hiện đại bằng dữ liệu genome toàn phần là một minh chứng điển hình cho sức mạnh của phương pháp này.

Phần mềm IQ-TREE là công cụ xây dựng cây tiến hoá dựa trên phương pháp **cực đại khả tín** (maximum likelihood), cho phép chọn mô hình tiến hoá tốt nhất và đánh giá độ tin cậy qua bootstrap. BEAST sử dụng **suy luận Bayes** (Bayesian inference) để xây dựng cây tiến hoá có thời gian, tích hợp đồng hồ phân tử và dữ liệu hoá thạch để ước tính thời điểm phân kỳ một cách thống kê.

Một thách thức quan trọng trong phylogenomics là hiện tượng **gene tree discordance** (bất đồng thuận cây gene: khi các gene khác nhau cho ra những cây tiến hoá không nhất quán do phân loại dòng (incomplete lineage sorting) hoặc dòng gene lịch sử. Các phương pháp liên kết loài (**species tree methods**) như ASTRAL được phát triển để tổng hợp các cây gene thành một cây loài thống nhất bằng cách mô hình hoá sự bất đồng thuận này.

### 3.3. DNA Cổ Đại

**DNA cổ đại** (ancient DNA, aDNA) là DNA được khai thác từ mẫu vật sinh học đã được bảo quản hàng chục đến hàng trăm nghìn năm, bao gồm xương, răng và mô đông lạnh từ lớp băng vĩnh cửu. Khả năng giải trình tự aDNA đã mở ra một cuộc cách mạng trong khảo cổ học và nhân học tiến hoá.

Thách thức kỹ thuật chính của aDNA là sự suy thoái của DNA theo thời gian: các phân tử ngắn đi, bị phá vỡ và bị biến đổi hoá học theo những cách đặc trưng (đặc biệt là sự chuyển đổi C thành T ở đầu phân tử). Phần mềm MapDamage và Mapdamage2 được dùng để nhận diện và hiệu chỉnh những dấu hiệu tổn thương này.

Nghiên cứu của Haak et al. (2015) về DNA cổ đại của hàng trăm cá thể từ thảo nguyên Á-Âu đã chứng minh rằng làn sóng di cư từ thảo nguyên Pontic-Caspian là nguồn gốc của cấu trúc gene châu Âu hiện đại và có thể liên quan đến sự phát tán của ngữ hệ Ấn-Âu, giải quyết một câu hỏi lịch sử kéo dài hàng thế kỷ.

### 3.4. Pangenomics

**Pangenomics** là cách tiếp cận đại diện cho genome của một loài không phải bằng một trình tự tham chiếu duy nhất, mà bằng một **đồ thị genome** (genome graph) bao trùm toàn bộ sự đa dạng di truyền trong loài đó. Khái niệm pangenome lần đầu được Tettelin et al. (2005) đề xuất cho vi khuẩn, phân biệt giữa **core genome** (core genome: tập hợp gene có mặt ở tất cả các chủng) và **accessory genome** (accessory genome: tập hợp gene chỉ có ở một số chủng).

Ngày nay, pangenomics đã được áp dụng cho thực vật, động vật và người. Pangenome người, được công bố năm 2023 bởi Human Pangenome Reference Consortium, tích hợp genome được lắp ráp chất lượng cao của 47 cá thể từ nhiều nhóm dân số khác nhau trên thế giới, cung cấp một bức tranh đầy đủ hơn về sự đa dạng cấu trúc genome người mà trình tự tham chiếu đơn lẻ không thể nắm bắt.

## 4. Ứng Dụng Thực Tiễn

### 4.1. Y Học và Sức Khoẻ Con Người

Hiểu lịch sử tiến hoá của genome người giúp làm sáng tỏ cơ chế của nhiều bệnh lý phức tạp. Chọn lọc tự nhiên đã định hình các immune gene, metabolic gene và neurological gene theo những cách có thể vừa bảo vệ vừa gây hại trong bối cảnh y tế hiện đại.

**Giả thuyết tiết kiệm** (thrifty gene hypothesis) đề xuất rằng các biến thể gene liên quan đến tích trữ năng lượng hiệu quả (bao gồm kháng insulin) được chọn lọc tích cực trong môi trường khan hiếm thức ăn hàng nghìn năm trước, nhưng lại trở thành yếu tố nguy cơ cho bệnh tiểu đường tuýp 2 trong điều kiện thừa dinh dưỡng ngày nay. Đây là ví dụ điển hình về **xung đột tiến hoá và y học** (evolutionary medicine).

Genomics tiến hoá còn giúp nghiên cứu sự tiến hoá của các mầm bệnh. Phân tích genome của virus, vi khuẩn và ký sinh trùng theo thời gian thực cho phép theo dõi sự xuất hiện và lan rộng của các biến chủng kháng thuốc, đặc biệt có giá trị trong kiểm soát dịch bệnh như đã thấy rõ trong đại dịch COVID-19 với hệ thống giám sát SARS-CoV-2 toàn cầu.

### 4.2. Bảo Tồn Đa Dạng Sinh Học

Genomics tiến hoá cung cấp những công cụ thiết yếu cho sinh học bảo tồn. Phân tích đa dạng di truyền trong các quần thể bị đe doạ giúp xác định mức độ suy giảm vốn gene, phát hiện các **đơn vị tiến hoá có tầm quan trọng bảo tồn** (evolutionarily significant units, ESUs), và hỗ trợ thiết kế chương trình nhân giống bảo tồn.

Nghiên cứu genome của các loài nguy cấp như tê giác Sumatra, báo tuyết, và cá heo sông Dương Tử đã sử dụng các phương pháp này để ước tính kích thước quần thể hiệu quả trong lịch sử, xác định mức độ cận huyết hiện tại và đánh giá khả năng thích nghi với biến đổi khí hậu. Kết quả trực tiếp hỗ trợ quyết định về chiến lược bảo tồn và ưu tiên phân bổ nguồn lực.

### 4.3. Nông Nghiệp và Chọn Giống

Hiểu nguồn gốc tiến hoá và quá trình thuần hoá của cây trồng vật nuôi giúp cải thiện chương trình chọn giống hiện đại. Quá trình thuần hoá để lại dấu ấn đặc trưng trong genome: **dấu ấn thuần hoá** (domestication sweep) là những vùng genome có đa dạng di truyền giảm mạnh do chọn lọc nhân tạo tập trung.

Phân tích genome lúa, ngô và lúa mì đã xác định những genomic region chịu chọn lọc trong quá trình thuần hoá, từ đó mở ra hướng cải tạo giống nhắm vào các đặc tính năng suất, chất lượng dinh dưỡng và khả năng chịu đựng stress môi trường. Trong chăn nuôi, genomics tiến hoá giúp truy tìm lịch sử phả hệ của các giống và đánh giá nguy cơ suy thoái do lai cận huyết, đặc biệt quan trọng với các giống địa phương có số lượng ít.

## Kết Luận

Genomics tiến hoá đã thay đổi cách chúng ta nhìn nhận sự sống. Bằng cách kết hợp sức mạnh của giải trình tự genome quy mô lớn với lý thuyết tiến hoá được kiểm chứng qua hơn 150 năm, lĩnh vực này đã trả lời nhiều câu hỏi lớn về nguồn gốc loài người, cơ chế thích nghi và lịch sử sự sống trên Trái Đất. Từ việc xác định gene kháng bệnh được chọn lọc qua hàng nghìn năm, đến việc tái dựng làn sóng di cư của người tiền sử, genomics tiến hoá kết nối quá khứ sâu với y học và sinh thái học của hiện tại.

Với sự phát triển của giải trình tự đọc dài (**long-read sequencing**) cho phép lắp ráp genome chất lượng cao, pangenomics đại diện cho toàn bộ đa dạng di truyền trong loài, và các phương pháp suy luận dựa trên học máy cho các bộ dữ liệu genomics quy mô lớn, lĩnh vực này đang tiếp tục mở rộng ranh giới của những gì chúng ta có thể biết về genome sinh vật. Trong các bài viết tiếp theo, chúng tôi sẽ đề cập quy trình thực hành phân tích genomics quần thể và xây dựng cây phát sinh loài từ dữ liệu thực.

## Tài Liệu Tham Khảo

Nielsen, R., Hellmann, I., Hubisz, M., Bustamante, C., & Clark, A. G. (2007). Recent and ongoing selection in the human genome.
*Nature Reviews Genetics*, *8*(11), 857–868. https://doi.org/10.1038/nrg2187

Jarvis, E. D., Mirarab, S., Aberer, A. J., Li, B., Houde, P., Li, C., Ho, S. Y. W., Faircloth, B. C., Nabholz, B., Howard, J. T., Suh, A., Weber, C. C., da Fonseca, R. R., Li, J., Zhang, F., Li, H., Zhou, L., Narula, N., Liu, L., ... Zhang, G. (2014). Whole-genome analyses resolve early branches in the tree of life of modern birds.
*Science*, *346*(6215), 1320–1331. https://doi.org/10.1126/science.1253451

Haak, W., Lazaridis, I., Patterson, N., Rohland, N., Mallick, S., Llamas, B., Brandt, G., Nordenfelt, S., Harney, E., Stewardson, K., Fu, Q., Mittnik, A., Bánffy, E., Economou, C., Francken, M., Friederich, S., Pena, R. G., Hallgren, F., Khartanovich, V., ... Reich, D. (2015). Massive migration from the steppe was a source for Indo-European languages in Europe.
*Nature*, *522*(7555), 207–211. https://doi.org/10.1038/nature14317

Tettelin, H., Masignani, V., Cieslewicz, M. J., Donati, C., Medini, D., Ward, N. L., Angiuoli, S. V., Crabtree, J., Jones, A. L., Durkin, A. S., DeBoy, R. T., Davidsen, T. M., Mora, M., Scarselli, M., Margarit, I., Peterson, J. D., Hauser, C. R., Sundaram, J. P., Nelson, W. C., ... Fraser, C. M. (2005). Genome analysis of multiple pathogenic isolates of Streptococcus agalactiae: implications for the microbial "pan-genome".
*Proceedings of the National Academy of Sciences*, *102*(39), 13950–13955. https://doi.org/10.1073/pnas.0506758102

Wohns, A. W., Wong, Y., Jeffery, B., Akbari, A., Mallick, S., Pinhasi, R., Patterson, N., Reich, D., Kelleher, J., & McVean, G. (2022). A unified genealogy of modern and ancient genomes.
*Science*, *375*(6583), eabi8264. https://doi.org/10.1126/science.abi8264

Human Pangenome Reference Consortium. (2023). A draft human pangenome reference.
*Nature*, *617*(7960), 312–324. https://doi.org/10.1038/s41586-023-05896-x
