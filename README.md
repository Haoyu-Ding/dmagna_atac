# dmagna_atac

独立的 *Daphnia magna* paired-end ATAC-seq Nextflow pipeline。第一版采用样本表作为唯一输入选择来源，不会扫描 FASTQ 目录，也不会根据文件名猜测样本身份。

## 当前实现

- 严格校验并规范化 CSV sample sheet；校验通过后才进入分析流程。
- 允许同一个 `sample_id` 对应多个 lane，并在 trimming 前合并对应 R1/R2。
- 原始 FASTQ FastQC。
- Cutadapt adapter/quality trimming。
- BWA-MEM2 paired-end alignment、排序和索引。
- `samtools fixmate`、`markdup`、proper-pair/低质量 read 过滤和线粒体 contig 去除。
- MACS3 `BAMPE` 模式逐样本 peak calling。
- trimmed-read QC、MultiQC 和 Nextflow execution report/timeline/trace。
- 记录 `samplesheet.used.csv`、`input_files.tsv`、输入校验报告和软件版本。

参考 FASTA、BWA-MEM2 index、线粒体 contig 名称和 MACS3 effective genome size 都必须由使用者针对最终采用的 *D. magna* reference 明确提供，因此没有写死在代码中。

## 先运行样本表校验测试

```bash
python3 -m unittest discover -s tests -v
```

## 准备容器

在可以构建 Docker 镜像的环境中：

```bash
docker build -t dmagna-atac:0.1.0 .
```

VSC 使用 Apptainer。可以把同一个镜像发布到一个 VSC 可访问的 registry，或者在有 Docker 的机器上转换成 `.sif` 后上传。`nextflow.config` 已提供 `docker` 和 `apptainer` 两个 profile。

## 准备参考基因组

需要先确定并记录最终采用的 *D. magna* reference 版本，随后准备 FASTA 和 BWA-MEM2 index。示意命令如下：

```bash
bwa-mem2 index -p /path/to/reference/bwa_index/dmagna.fa /path/to/reference/dmagna.fa
```

`--bwa_index_dir` 应指向包含该 prefix 产生的全部 index 文件的目录。线粒体 contig 名称和 `--macs_genome_size` 必须来自同一个 reference 配置，不能凭样本文件名推断。

## VSC 运行命令

在 VSC 上 clone 仓库、准备 sample sheet 和 reference 后：

```bash
nextflow run . \
  -profile apptainer \
  -resume \
  --samplesheet /path/to/samplesheet.csv \
  --reference_fasta /path/to/reference/dmagna.fa \
  --bwa_index_dir /path/to/reference/bwa_index \
  --mitochondrial_contig '<MITO_CONTIG_FROM_REFERENCE>' \
  --macs_genome_size '<EFFECTIVE_GENOME_SIZE>' \
  --outdir /path/to/results
```

第一次正式运行前，建议先用一个小的、真实的 paired-end 子集验证容器、参考索引、线粒体 contig 和 scheduler 配置。正式运行应保留 `results/pipeline_info/` 中的输入清单和执行记录。

## 输出目录

主要输出位于：

```text
results/
├── qc/raw/<sample_id>/
├── qc/trimmed/<sample_id>/
├── qc/multiqc/
├── alignment/<sample_id>.sorted.bam[.bai]
├── filtered_bam/<sample_id>.filtered.bam[.bai]
├── peaks/<sample_id>.narrowPeak
└── pipeline_info/
```

## 下一阶段

当前代码已经覆盖第一版的工程骨架，但正式分析前仍需完成三件事：

1. 从实验记录确认四个样本的 `condition`、`replicate`，并把真实 FASTQ 路径写入 sample sheet。
2. 确认 reference 版本、线粒体 contig 名称和 MACS3 effective genome size，并在 VSC 生成/验证 BWA-MEM2 index。
3. 在 VSC 先跑小规模 smoke test，再跑四样本完整流程；检查 `input_validation_report.txt`、MultiQC、BAM 统计和每个样本的 peak 文件。

差异可及性、共识 peak、peak 注释、motif 和功能富集不属于当前版本，等四个样本的 QC 和生物学重复关系确认后再单独设计。
