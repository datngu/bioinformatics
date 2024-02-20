---
layout: post
title:  "Giới thiệu về Terminal, Bash Script và những dòng lệnh cơ bản"
author: dat
categories: [ Data science, Bioinformatics ]
image: assets/my_figs/ds/bash_script.jpg

---

Terminal/ Bash Script được sử dụng cực kỳ phổ biến trong tin sinh học, chủ yếu dùng để thao tác và điều khiển những máy tính có hiệu năng cao (HPC) để thực hiện những tính toán chuyên biệt. Trong bài viết này, chúng tôi sẽ giới thiệu sơ lược về Terminal và Bash Script cho bạn nào mới/ muốn bắt đầu với tin sinh học!



**Terminal:**

* Terminal (còn gọi là Command Line Interface - CLI) là giao diện dòng lệnh cho phép người dùng tương tác trực tiếp với hệ thống bằng cách nhập các lệnh văn bản.
* Terminal cung cấp khả năng truy cập và điều khiển hệ thống mạnh mẽ, linh hoạt hơn so với giao diện đồ họa (GUI).
* Terminal thường được sử dụng bởi các quản trị viên hệ thống, lập trình viên và những người dùng có kinh nghiệm về hệ thống.
* Ngoài ra, Terminal cũng được sử dụng cực kỳ phổ biến trong tin sinh học, chủ yếu dùng để thao tác và điều khiển những máy tính có hiệu năng cao (HPC) để thực hiện những tính toán chuyên biệt.

**Bash Script:**

Thay vì phải copy và paste từng dòng lệnh vào terminal, và lặp đi lặp lại, Bash Scipt là bước tiếp theo để giúp bạn liệt kê những câu lệnh vào 1 file text, và giúp máy tính thực hiện công việc một cách tự động...

* Bash Script là tập tin văn bản chứa các lệnh bash được thực thi tuần tự. 
* Bash Script giúp tự động hóa các tác vụ lặp đi lặp lại trên hệ thống, tiết kiệm thời gian và công sức.
* Bash Script có thể được sử dụng để thực hiện nhiều chức năng khác nhau như:
    * Tự động hóa các tác vụ quản trị hệ thống
    * Xử lý dữ liệu
    * Viết chương trình đơn giản
    * Tạo các công cụ tùy chỉnh

**Lợi ích của việc sử dụng Bash Script:**

* **Tự động hóa:** Giúp thực hiện các tác vụ lặp đi lặp lại một cách tự động, tiết kiệm thời gian và công sức.
* **Nâng cao hiệu quả:** Cho phép thực hiện nhiều tác vụ cùng lúc, tăng hiệu quả công việc.
* **Giảm thiểu lỗi:** Giúp giảm thiểu lỗi do thao tác thủ công.
* **Tăng khả năng tùy chỉnh:** Cho phép tùy chỉnh các tác vụ theo nhu cầu cụ thể.

**Cấu trúc cơ bản của Bash Script:**

* **Dòng Shebang:** Dòng đầu tiên của tập tin bắt đầu với `#!` và chỉ định trình thông dịch shell sẽ sử dụng để thực thi tập tin. Ví dụ: `#!/bin/bash` sẽ chạy chương trình bằng `bash`, `#!/bin/python3` sẽ chạy chương trình bằng `python3`.
* **Bình luận:** Dòng bắt đầu với `#` là bình luận và không được thực thi.
* **Lệnh:** Các dòng tiếp theo chứa các lệnh bash sẽ được thực thi tuần tự.
* **Biến:** Biến được sử dụng để lưu trữ giá trị và có thể được sử dụng trong các lệnh.
* **Điều kiện:** Câu lệnh điều kiện cho phép thực thi các lệnh khác nhau dựa trên điều kiện cụ thể.
* **Vòng lặp:** Vòng lặp cho phép thực thi một khối lệnh nhiều lần.

**Một số câu lệnh bash cơ bản:**

* **ls:** Liệt kê các tập tin trong thư mục hiện tại.
* **cd:** Thay đổi thư mục hiện tại.
* **pwd:** Hiển thị đường dẫn đến thư mục hiện tại.
* **mkdir:** Tạo thư mục mới.
* **rmdir:** Xóa thư mục rỗng.
* **cp:** Sao chép tập tin.
* **mv:** Di chuyển tập tin.
* **rm:** Xóa tập tin.
* **echo:** Hiển thị thông tin ra màn hình.
* **cat:** Hiển thị nội dung tập tin.
* **less:** Hiển thị nội dung tập tin nhưng không in ra terminal.
* **grep:** Tìm kiếm chuỗi văn bản trong tập tin.
* **sort:** Sắp xếp các dòng trong tập tin.
* **uniq:** Loại bỏ các dòng trùng lặp trong tập tin.
* **head:** Hiển thị các dòng đầu tiên của tập tin.
* **tail:** Hiển thị các dòng cuối cùng của tập tin.


**Ví dụ về Bash Script đơn giản:**

```bash
#!/bin/bash

# Tạo thư mục mới
mkdir new_folder

# Di chuyển vào thư mục mới
cd new_folder

# Tạo tập tin mới
echo "Hello World!" > hello_world.txt

# Hiển thị nội dung tập tin
cat hello_world.txt

# Xóa thư mục mới
rmdir new_folder
```

**Lưu ý:**

* Khi viết Bash Script, cần chú ý đến cú pháp và định dạng.
* Nên sử dụng bình luận để giải thích các lệnh trong tập tin.
* Nên kiểm tra kỹ lưỡng tập tin trước khi thực thi.


**Để tìm hiểu thêm về Bash Script:**

* Tham khảo tài liệu trực tuyến:
    * [https://www.tldp.org/LDP/abs/html/](https://www.tldp.org/LDP/abs/html/)


**Kết luận:**

* Terminal, Bash Script và các dòng lệnh cơ bản là những công cụ quan trọng giúp bạn thao tác và quản trị hệ thống hiệu quả.
* Nắm vững những công cụ này giúp bạn tăng năng suất và khả năng kiểm soát hệ thống cũng như thực hiện phân tích tin sinh học hiệu quả hơn.



