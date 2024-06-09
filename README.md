# Final-Project-Construction-and-implementation-of-data-analysis-process
Group Members: 赵筱语 仝晓亮

[toc]

## 1. 摘要：

下载果蝇gtf注释文件、hisat2索引文件和Dif1突变果蝇在孵化后108 小时和132小时的大脑的RNA-seq数据，然后对RNA-seq数据进行数据质量控制、质量过滤和去除接头、数据比对、表达定量、差异表达分析和数据可视化，从而直观了解Dif1 突变果蝇在特定时间点的基因表达情况。



## 2. 前言：

1）背景介绍：Dif1是果蝇的一个关键基因，编码一种转录因子，参与调控果蝇的免疫应答和抗菌活动。研究Dif1突变果蝇可以帮助我们深入了解基因在免疫调控、发育和其他生物学过程中的作用机制。这有助于深入了解免疫相关疾病的发病机理，并为新药研发提供线索。

2）工作思路：首先下载RNA-seq分析所需工具，然后从NCBI查找不同条件下的RNA-seq数据编号，每个条件下载两个sra文件(prefetch)，并下载参考基因组和索引文件，对不同条件下的RNA-seq数据进行数据质量控制(fastqc)、质量过滤和去除接头(Trimmomatic)、数据比对(HISAT2)、表达定量(featureCounts)、差异表达分析(DESeq2)和数据可视化（热图、MA图、火山图）。

3）工作概要：

① 数据准备：

下载果蝇的gtf注释文件、hisat2索引文件和Dif1突变果蝇在不同时间点的RNA-seq数据。

② 数据处理：

进行数据质量控制和质量过滤，包括使用FastQC评估数据质量和使用Trimmomatic去除接头和低质量序列。

③ 数据比对：

使用Hisat2或Bowtie2进行RNA-seq数据与参考基因组的比对，生成比对结果。

④ 表达定量：

使用featureCounts或HTSeq对比对后的数据进行表达定量，得到基因的计数信息。

⑤ 差异表达分析：

使用DESeq2、edgeR或limma等工具进行差异表达分析，比较不同时间点的基因表达情况，找出差异表达的基因。

⑥ 数据可视化：

使用R包（如ggplot2）绘制热图、Volcano图等，可视化差异表达基因的结果，从而直观了解基因表达情况。

⑦ 结果解释：

综合分析差异表达基因的结果，找出在不同时间点上调和下调的基因，探索Dif1突变果蝇在基因表达水平上的差异，以及在特定时间点的生物学特点。



## 3. 数据集与方法：

1）项目涉及的工具：

① samtools：用prefetch从NCBI Sequence Read Archive（SRA）数据库中下载SRA格式的测序数据文件。

② SRA Toolkit：用fastq-dump将SRA文件转换成FASTQ格式的测序数据文件，并可以选择按照read pairs进行拆分（--split-files选项）。

③ FastQC：用fastqc对FASTQ文件进行测序数据的质量控制分析。

④ Trimmomatic：一个用于对测序数据进行质量过滤和去除接头的工具。

⑤ HISAT2：用于进行RNA-Seq数据的比对。

⑥ featureCounts：用于进行基因表达定量分析，根据比对结果统计每个基因的reads数量。

⑦ DBI：R语言中用于数据库操作的库。

⑧ rlang：R语言中用于语言解析和操作的库。

⑨ DESeq2：用于差异表达分析的R包。

⑩ pheatmap：用于绘制热图的R包。

⑪ ggplot2：用于绘制数据可视化图形的R包。

2）项目涉及的数据：

① Dif1 突变果蝇在孵化后 108 小时的大脑的RNA-seq数据：SRR29185161，SRR29185162

② Dif1 突变果蝇在孵化后 132 小时的大脑的RNA-seq数据：SRR29185159，SRR29185160

③ gtf注释文件：Drosophila_melanogaster.BDGP6.32.107.chr.gtf

④ hisat2索引文件：bdgp6.tar.gz（解压后改名为hisat2-index）

3）相关工具下载方法：

① samtools：

```
$ mkdir $HOME/bin
$ export PATH=$HOME/bin:$PATH
$ wget https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2
$ tar jxvf samtools-1.20.tar.bz2 
$ cd samtools-1.20
$ ./configure
$ make
$ sudo make install
$ cd
$ cp samtools-1.20/samtools $HOME/bin
```

② SRA Toolkit：

```
$ wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
$ tar -xzvf sra-toolkit.tar.gz
$ cd sratoolkit.current-ubuntu64
$ ./installer.sh
```

③ FastQC：

```swift
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
unzip fastqc_v0.11.9.zip
cd FastQC
chmod +x fastqc
./fastqc
```

④ Trimmomatic：

```
$ wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
$ unzip Trimmomatic-0.39.zip
```

⑤ HISAT2：

```
$ sudo apt install hisat2
```

⑥ featureCounts：

```
$ sudo apt install subread
```

⑦ DBI：

```
> install.packages("DBI")
> library(DBI)
```

⑧ rlang：

```
> install.packages("rlang")
> library(rlang)
```

⑨ DESeq2：

```
> if (!requireNamespace("BiocManager", quietly = TRUE))
+ install.pack.packages("BiocManager")
> BiocManager::install("DESeq2")
> library(DESeq2)
```

⑩ pheatmap：

```
> install.packages("pheatmap")
> library(pheatmap)
```

⑪ ggplot2：

```
> install.packages("ggplot2")
> library(ggplot2)
```

4）详细的分析流程：

① 下载注释文件和索引文件：

```
# 创建目录存放分析数据：
$ mkdir my_rnaseq_exp/
$ cd my_rnaseq_exp/

# 下载gtf注释文件：
$ wget http://ftp.ensembl.org/pub/release-107/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.32.107.chr.gtf.gz

# 下载RNA-Seq需要的hisat2的索引文件：
$ wget https://genome-idx.s3.amazonaws.com/hisat/bdgp6.tar.gz
# 解压下载的索引文件：
$ tar -zxvf bdgp6.tar.gz
# 改名：
$ mv bdgp6 hisat2-index
```

② 下载RNA-seq数据并转为fastq格式：

```
# 下载Dif1 突变果蝇在孵化 132 小时后的大脑的RNA-seq数据：
$ prefetch SRR29185159.sra SRR29185160.sra
# 下载Dif1 突变果蝇在孵化 108 小时后的大脑的RNA-seq数据：
$ prefetch SRR29185161.sra SRR29185162.sra

创建目录存放fastq格式的RNA-seq数据：
$ mkdir raw_fq1/
$ mkdir raw_fq2/
$ mkdir raw_fq3/
$ mkdir raw_fq4/

# 将sra格式的RNA-seq数据转为fastq格式：
$ fastq-dump --split-files -O /home/zxy/my_rnaseq_exp/raw_fq1 /home/zxy/my_rnaseq_exp/SRR29185159/SRR29185159.sra
$ fastq-dump --split-files -O /home/zxy/my_rnaseq_exp/raw_fq2 /home/zxy/my_rnaseq_exp/SRR29185160/SRR29185160.sra
$ fastq-dump --split-files -O /home/zxy/my_rnaseq_exp/raw_fq3 /home/zxy/my_rnaseq_exp/SRR29185161/SRR29185161.sra
$ fastq-dump --split-files -O /home/zxy/my_rnaseq_exp/raw_fq4 /home/zxy/my_rnaseq_exp/SRR29185162/SRR29185162.sra
```

③ 数据质量控制和过滤：

```
# 用fastqc对数据进行数据质量控制：
$ fastqc /home/zxy/my_rnaseq_exp/raw_fq1/SRR29185159_1.fastq
$ fastqc /home/zxy/my_rnaseq_exp/raw_fq1/SRR29185159_2.fastq

$ fastqc /home/zxy/my_rnaseq_exp/raw_fq2/SRR29185160_1.fastq
$ fastqc /home/zxy/my_rnaseq_exp/raw_fq2/SRR29185160_2.fastq

$ fastqc /home/zxy/my_rnaseq_exp/raw_fq3/SRR29185161_1.fastq
$ fastqc /home/zxy/my_rnaseq_exp/raw_fq3/SRR29185161_2.fastq

$ fastqc /home/zxy/my_rnaseq_exp/raw_fq4/SRR29185162_1.fastq
$ fastqc /home/zxy/my_rnaseq_exp/raw_fq4/SRR29185162_2.fastq

# 使用Trimmomatic对数据进行质量过滤和去除接头：
$ java -jar Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 /home/zxy/my_rnaseq_exp/raw_fq1/SRR29185159_1.fastq /home/zxy/my_rnaseq_exp/raw_fq1/SRR29185159_2.fastq /home/zxy/my_rnaseq_exp/raw_fq1/SRR29185159_1_paired.fastq /home/zxy/my_rnaseq_exp/raw_fq1/SRR29185159_1_unpaired.fastq /home/zxy/my_rnaseq_exp/raw_fq1/SRR29185159_2_paired.fastq /home/zxy/my_rnaseq_exp/raw_fq1/SRR29185159_2_unpaired.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

$ java -jar Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 /home/zxy/my_rnaseq_exp/raw_fq2/SRR29185160_1.fastq /home/zxy/my_rnaseq_exp/raw_fq2/SRR29185160_2.fastq /home/zxy/my_rnaseq_exp/raw_fq2/SRR29185160_1_paired.fastq /home/zxy/my_rnaseq_exp/raw_fq2/SRR29185160_1_unpaired.fastq /home/zxy/my_rnaseq_exp/raw_fq2/SRR29185160_2_paired.fastq /home/zxy/my_rnaseq_exp/raw_fq2/SRR29185160_2_unpaired.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

$ java -jar Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 /home/zxy/my_rnaseq_exp/raw_fq3/SRR29185161_1.fastq /home/zxy/my_rnaseq_exp/raw_fq3/SRR29185161_2.fastq /home/zxy/my_rnaseq_exp/raw_fq3/SRR29185161_1_paired.fastq /home/zxy/my_rnaseq_exp/raw_fq3/SRR29185161_1_unpaired.fastq /home/zxy/my_rnaseq_exp/raw_fq3/SRR29185161_2_paired.fastq /home/zxy/my_rnaseq_exp/raw_fq3/SRR29185161_2_unpaired.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

$ java -jar Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 /home/zxy/my_rnaseq_exp/raw_fq4/SRR29185162_1.fastq /home/zxy/my_rnaseq_exp/raw_fq4/SRR29185162_2.fastq /home/zxy/my_rnaseq_exp/raw_fq4/SRR29185162_1_paired.fastq /home/zxy/my_rnaseq_exp/raw_fq4/SRR29185162_1_unpaired.fastq /home/zxy/my_rnaseq_exp/raw_fq4/SRR29185162_2_paired.fastq /home/zxy/my_rnaseq_exp/raw_fq4/SRR29185162_2_unpaired.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

④ 数据比对：

```
# 使用HISAT2进行数据比对：
$ hisat2 -x /home/zxy/my_rnaseq_exp/hisat2-index/genome -1 /home/zxy/my_rnaseq_exp/raw_fq1/SRR29185159_1_paired.fastq -2 /home/zxy/my_rnaseq_exp/raw_fq1/SRR29185159_2_paired.fastq -S /home/zxy/my_rnaseq_exp/SRR29185159.sam

$ hisat2 -x /home/zxy/my_rnaseq_exp/hisat2-index/genome -1 /home/zxy/my_rnaseq_exp/raw_fq2/SRR29185160_1_paired.fastq -2 /home/zxy/my_rnaseq_exp/raw_fq2/SRR29185160_2_paired.fastq -S /home/zxy/my_rnaseq_exp/SRR29185160.sam

$ hisat2 -x /home/zxy/my_rnaseq_exp/hisat2-index/genome -1 /home/zxy/my_rnaseq_exp/raw_fq3/SRR29185161_1_paired.fastq -2 /home/zxy/my_rnaseq_exp/raw_fq3/SRR29185161_2_paired.fastq -S /home/zxy/my_rnaseq_exp/SRR29185161.sam

$ hisat2 -x /home/zxy/my_rnaseq_exp/hisat2-index/genome -1 /home/zxy/my_rnaseq_exp/raw_fq4/SRR29185162_1_paired.fastq -2 /home/zxy/my_rnaseq_exp/raw_fq4/SRR29185162_2_paired.fastq -S /home/zxy/my_rnaseq_exp/SRR29185162.sam
```

⑥ 表达定量：

```
用featureCounts进行表达定量：
$ featureCounts -a /home/zxy/my_rnaseq_exp/Drosophila_melanogaster.BDGP6.32.107.chr.gtf -o /home/zxy/my_rnaseq_exp/SRR29185159_counts.txt -p /home/zxy/my_rnaseq_exp/SRR29185159.sam

$ featureCounts -a /home/zxy/my_rnaseq_exp/Drosophila_melanogaster.BDGP6.32.107.chr.gtf -o /home/zxy/my_rnaseq_exp/SRR29185160_counts.txt -p /home/zxy/my_rnaseq_exp/SRR29185160.sam

$ featureCounts -a /home/zxy/my_rnaseq_exp/Drosophila_melanogaster.BDGP6.32.107.chr.gtf -o /home/zxy/my_rnaseq_exp/SRR29185161_counts.txt -p /home/zxy/my_rnaseq_exp/SRR29185161.sam

$ featureCounts -a /home/zxy/my_rnaseq_exp/Drosophila_melanogaster.BDGP6.32.107.chr.gtf -o /home/zxy/my_rnaseq_exp/SRR29185162_counts.txt -p /home/zxy/my_rnaseq_exp/SRR29185162.sam
```

⑦ 差异表达分析：

```
$ sudo R

> library(DBI)
> library(rlang)
> library(DESeq2)

# 132h
> SRR29185159 <- as.matrix(read.table("SRR29185159_counts.txt", header = TRUE, row.names = 1))
> SRR29185160 <- as.matrix(read.table("SRR29185160_counts.txt", header = TRUE, row.names = 1))

> countData1 <- as.matrix(SRR29185159[, ncol(SRR29185159)])
> countData2 <- as.matrix(SRR29185160[, ncol(SRR29185160)])
> countData1 <- apply(countData1, 2, as.numeric)
> countData2 <- apply(countData2, 2, as.numeric)

> countData_132 <- cbind(countData1, countData2)
> colnames(countData_132) <- c("SRR29185159", "SRR29185160")

# 108h
> SRR29185161 <- as.matrix(read.table("SRR29185161_counts.txt", header = TRUE, row.names = 1))
> SRR29185162 <- as.matrix(read.table("SRR29185162_counts.txt", header = TRUE, row.names = 1))

> countData3 <- as.matrix(SRR29185161[, ncol(SRR29185161)])
> countData4 <- as.matrix(SRR29185162[, ncol(SRR29185162)])
> countData3 <- apply(countData3, 2, as.numeric)
> countData4 <- apply(countData4, 2, as.numeric)

> countData_108 <- cbind(countData3, countData4)
> colnames(countData_108) <- c("SRR29185161", "SRR29185162")

# 合并数据
> countData <- cbind(countData_132, countData_108)
> rownames(countData) <- rownames(SRR29185159)

# 把合并数据countData写到CSV文件里
write.csv(countData, file = "merged_count_data.csv")

# 创建样本信息数据框
> sampleInfo <- data.frame(
  row.names = colnames(countData),
  condition = factor(c(rep("132h", ncol(countData_132)), rep("108h", ncol(countData_108))))
)

# 创建 DESeqDataSet 对象
> dds <- DESeqDataSetFromMatrix(countData = countData, colData = sampleInfo, design = ~ condition)

# 进行差异表达分析
> dds <- DESeq(dds)
> res <- results(dds)
```

⑧ 数据可视化：

```
> library(pheatmap)
> library(ggplot2)

# 提取差异表达基因
> resOrdered <- res[order(res$padj),]

# 选择前100个差异表达基因用于绘制热图
> topGenes <- rownames(resOrdered)[1:100]

# 提取这些基因在每个样本中的表达值
> mat <- assay(rlog(dds))[topGenes,]

# 绘制热图
> pheatmap(mat, cluster_rows = TRUE, cluster_cols = TRUE)

# 绘制MA图
> plotMA(res, ylim = c(-2, 2))

# 去除NA值
> res <- na.omit(res)

# 绘制火山图
> alpha <- 0.01
> sig <- res[res$padj < alpha, ]
> not_sig <- res[res$padj >= alpha, ]
> plot <- ggplot(res, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(data = not_sig, aes(color = "black"), alpha = 0.6) +
  geom_point(data = sig, aes(color = "red"), alpha = 0.6) +
  theme_minimal() +
  scale_color_manual(values = c("black", "red")) +
  labs(x = "Log2 Fold Change", y = "-Log10 P-value") +
  geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  theme(legend.position = "none")

> print(plot)
```



## 4. 结果：数据分析结果展示

1）火山图：

![image-20240607125026801](C:\Users\ZHAO\AppData\Roaming\Typora\typora-user-images\image-20240607125026801.png)

2）MA图：

![image-20240607125059391](C:\Users\ZHAO\AppData\Roaming\Typora\typora-user-images\image-20240607125059391.png)

3）热图：

![image-20240607125119002](C:\Users\ZHAO\AppData\Roaming\Typora\typora-user-images\image-20240607125119002.png)

4）由`merged_count_data.csv` 和热图可知，

① FBgn0000064、FBgn0004507等基因在132h的表达量相对108h有所下调

② FBgn0000639、FBgn0085276等基因在132h的表达量相对108h有所上调

通过分析基因表达量随时间的变化，有助于深入理解生物体内部的基因调控机制、生物学过程的时间动态性以及基因在不同时空下的功能和作用，为生物学研究提供重要的参考和指导。



## 5. 讨论：

1）关键点：

① 正确查找并下载所需RNA-seq数据（每个条件至少下载两组数据）、注释文件和索引文件

② 进行数据处理（质量控制、过滤、比对、表达定量、差异表达分析）

③ 用R包对差异表达基因进行可视化，以直观了解基因表达情况。

④ 综合分析差异表达基因的结果，找出Dif1突变果蝇在不同时间点上调和下调的基因。

2）不足：

① 仅有Dif1突变果蝇的RNA-seq数据，缺少野生型RNA-seq数据作对照

② 时间变量较少，仅有Dif1突变果蝇在孵化后108 小时和132小时的大脑的RNA-seq数据，缺少其他时间的数据，难以确定基因表达量随时间的变化情况

③ 每个条件下的RNA-seq文件数量少，应尽量增加文件数量，以减少偶然误差



## 6. 参考文献

1）Trapnell, C., Roberts, A., Goff, L. *et al.* Differential gene and transcript expression analysis of RNA-seq experiments with TopHat and Cufflinks. *Nat Protoc* **7**, 562–578 (2012). https://doi.org/10.1038/nprot.2012.016

2）Pertea, M., Kim, D., Pertea, G. *et al.* Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown. *Nat Protoc* **11**, 1650–1667 (2016). https://doi.org/10.1038/nprot.2016.095

3）https://blog.csdn.net/yangl7/article/details/108880855

4）【The pipeline of RNA-seq（菜鸟教学）】https://mbd.baidu.com/ma/s/wkbvhr59

5）https://www.jianshu.com/p/37596ccc5ae7

6）https://www.jianshu.com/p/12bb59611aa4

7）https://www.jianshu.com/p/5a64cb2acc82

8）https://blog.csdn.net/weixin_39900139/article/details/109005162

9）https://cloud.tencent.com/developer/article/2071040



## 7. 附录：

1）docker配置：

```
# docker极速安装，快速完成不浪费时间
curl  -sSfL get.docker.io -o get_docker.sh
bash  get_docker.sh --mirror Aliyun

# 安装docker-compose
# 从以下网址下载docker-compose，将docker-compose文件放在path变量目录下如：/usr/local/bin并增加可执行权限
sudo chmod +x /usr/local/bin/docker-compose
https://github.com/docker/compose/releases

# 从Ubuntu20.04镜像开始构建，Ubuntu22.04据说会不定期杀掉占用资源过多的进程
FROM        ubuntu:20.04
# 1.设置账户字符编码为C.UTF-8，提高兼容性；设置时区，ssh登录密码为20201110（账户为root），连接端口为9022
ENV         LANG=C.UTF-8 TZ=Asia/Shanghai PS=20201110 port=9022
# 设置时区、更新镜像软件、安装aria2（下载工具替代wget，curl以获取更快的下载速度，容错/下载会自动重试）
# openssh服务并更新配置文件，使root账户可以登录、更新root账户密码为设置值
RUN         ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone && \
            apt update && apt upgrade -y && apt install openssh-server aria2 -y && apt autoremove -y
RUN         sed -i "s/#PermitEmptyPasswords no/PermitEmptyPasswords no/g" /etc/ssh/sshd_config && \
            sed -i "s/#PermitRootLogin prohibit-password/PermitRootLogin yes/g" /etc/ssh/sshd_config && \
            sed -i "s/#Port 22/Port ${port}/g" /etc/ssh/sshd_config && \
            sed -i "s/#ListenAddress 0.0.0.0/ListenAddress 0.0.0.0/g" /etc/ssh/sshd_config && \
            sed -i "s/#LoginGraceTime 2m/LoginGraceTime 2m/g" /etc/ssh/sshd_config && \
            echo root:${PS} | chpasswd
# aria2下载miniconda、安装、添加channel、删除安装文件（减小镜像体积）、初始化miniconda
RUN         aria2c https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -d ~ && \
            bash ~/Miniconda3-latest-Linux-x86_64.sh -b && \
            ~/miniconda3/bin/conda config --add channels defaults    && \
            ~/miniconda3/bin/conda config --add channels bioconda    && \
            ~/miniconda3/bin/conda config --add channels conda-forge && \
            rm -f ~/Miniconda3-latest-Linux-x86_64.sh && \
            ~/miniconda3/bin/conda init
# 复制该文件到镜像root目录下，condarc为清华源配置文件，国内提速可以注销该行
# COPY       --chown=root:root ./condarc /root/.condarc
# 暴露ssh连接端口
EXPOSE      $port
# 初始化镜像运行：根据配置项变量PS修改root密码，该密码可以运行时重新设置初始化，最后启动ssh服务
ENTRYPOINT  ["/bin/bash", "-c", "echo root:${PS} | chpasswd && service ssh start -D"]

# 将condarc文件和dockerfile放在同一目录下，构建镜像
docker     build -t doujiangbaozi/sliverworkspace:latest .

# 或者使用已经构建好的镜像，直接拉取到本地
docker     pull     doujiangbaozi/sliverworkspace:latest

channels:

  - conda-forge
  - bioconda
  - defaults
    show_channel_urls: true
    default_channels:
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r
    custom_channels:
      conda-forge: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
      msys2: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
      bioconda: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
      menpo: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
      pytorch: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud

version: "3"
services:
  TumorOnly:
    image: doujiangbaozi/sliverworkspace:latest
    container_name: TumorOnly
    volumes:
      - /media/sliver/Data/data:/opt/data:rw                               #挂载输入数据目录
      - /media/sliver/Manufacture/tumor-only/envs:/root/miniconda3/envs:rw #挂载envs目录
      - /media/sliver/Manufacture/tumor-only/config:/opt/config:rw         #挂载config目录
      - /media/sliver/Manufacture/sliver/ref:/opt/ref:rw                   #挂载reference目录
      - /media/sliver/Manufacture/tumor-only/result:/opt/result:rw         #挂载中间文件和结果目录
    ports:
      - "9022:9022"                                                        #ssh连接端口，此处可以映射为其他端口
    environment:
      - TZ=Asia/Shanghai                                                   #设置时区
      - PS=20191124                                                        #设置ssh密码

# 查看docker容器运行状态
docker ps

# 或者docker-compose.yml目录下运行
docker-compose ps

# 用到的环境变量，以最简单的fastqc，multiqc为例
export env=/root/miniconda3/envs  #conda环境软件安装目录，最好挂载物理机volume
export conf=/opt/config           #conda环境配置文件目录
export sn=RD1703007FFP            #样本编号，sample  number
export pn=TumorOnly               #项目编号，project number
export result=/opt/result         #中间文件输出目录

# conda检测环境是否存在，首次运行不存在创建该环境并安装软件
if [ ! -d "${envs}/${pn}.fastqc" ]; then
  conda env create -f ${conf}/fastqc.yaml
fi
# 切换到环境下运行fastqc、multiqc
conda activate "${pn}.fastqc"
mkdir  -p ${result}/${sn}/qc
fastqc -o ${result}/${sn}/qc \
	${result}/${sn}/trimmed/${sn}_trimmed_R1.fastq.gz \
    ${result}/${sn}/trimmed/${sn}_trimmed_R2.fastq.gz
multiqc ${result}/${sn}/ -f -o ${result}/${sn}/qc
# 退出环境
conda deactivate

name: TumorOnly.fastqc
channels:

  - conda-forge
  - bioconda
  - defaults
    dependencies:
  - fastqc=0.11.9
  - multiqc=1.13a
```

2）github链接：





## 8. 成员分工：

1）查找文献：仝晓亮、赵筱语

2）工具下载、RNA-seq数据下载、数据处理、分析、绘图：赵筱语

3）docker配置：仝晓亮

4）项目报告：赵筱语 









