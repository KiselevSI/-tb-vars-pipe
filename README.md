
## Описание

Пайплайн предназначен для автоматизации биоинформатического анализа, включая этапы контроля качества, выравнивания, анализа покрытия и вызова генетических вариантов из FASTQ/FASTA файлов или данных SRA.

## Зависимости

Перед запуском скрипта убедитесь, что установлены следующие программы:
- `fastqc`
- `bwa`
- `samtools`
- `mosdepth`
- `multiqc`
- `bcftools`
- `wget`
- `fastq-dump`

## Установка через conda

Установите необходимые программы с помощью conda:

```bash
conda create -n bioinf_pipeline -c bioconda fastqc bwa samtools mosdepth multiqc bcftools sra-tools
conda activate bioinf_pipeline
```

## Использование скрипта

Сделайте скрипт исполняемым и запустите:

```bash
chmod +x pipeline.sh
./pipeline.sh -i <input_folder_or_files_or_txt> -o <output_directory> -t <concurrency>
```

### Параметры
- `-i` – Путь к FASTQ/FASTA файлам, текстовый файл со списком SRA-акцессий или файл со списком путей.
- `-o` – Директория для результатов.
- `-t` – Количество параллельных задач.

### Примеры запуска

```bash
# Папка с FASTQ-файлами
./pipeline.sh -i data/FASTQ/ -o results -t 4

# Файл со списком FASTQ-файлов
./pipeline.sh -i list_of_fastq.txt -o results -t 4

# Файл со списком SRA-акцессий
./pipeline.sh -i list_of_sra_accessions.txt -o results -t 4

# Отдельные FASTQ-файлы
./pipeline.sh -i sample1.fastq sample2.fastq -o results -t 4
```

## Структура результатов

```
<output_directory>
├── BAM
├── BCF
├── FASTQ
├── VCF
├── ref
└── stats
    ├── fastqc
    ├── mosdepth
    ├── multiqc
    └── samtools_stats
```

## Этапы пайплайна
- Скачивание данных из SRA (при необходимости)
- Контроль качества (FastQC)
- Выравнивание на референс (`bwa mem`)
- Анализ покрытия (mosdepth)
- Вызов вариантов (bcftools)
- Формирование итогового отчёта (multiqc)

## Итоговый отчёт

Итоговые отчёты находятся в папке:

```
<output_directory>/stats/multiqc
```



