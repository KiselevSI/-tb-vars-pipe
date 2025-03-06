#!/bin/bash
# Биоинформатический пайплайн
#
# Использование:
#   ./pipeline.sh -i <input_folder_or_files_or_txt> -o <output_directory> -t <concurrency>
#
# Параметры:
#   -i    Папка с FASTQ/FASTA файлами, текстовый файл со списком SRA-акцессий или список файлов (например, results12/FASTQ/*)
#   -o    Каталог для хранения всех результатов (будут созданы подпапки: FASTQ, ref, BAM, stats и т.д.)
#   -t    Количество одновременно обрабатываемых сэмплов (параллелизм на уровне задач)
#
# Параметры программ (количество потоков):
CORES_FASTQC=2
CORES_BWA=2
CORES_SAMTOOLS=2
CORES_MOSDEPTH=2
CORES_BCFT=2

#########################
# Проверка наличия основных команд
for cmd in fastqc bwa samtools mosdepth multiqc bcftools wget fastq-dump; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "Ошибка: $cmd не найден. Установите $cmd и убедитесь, что он доступен в PATH."
    exit 1
  fi
done

#########################
# Функция для скачивания данных из SRA.
download_sra() {
    local accession="$1"
    local outdir="$2"
    echo "Скачивание SRA: $accession..."
    fastq-dump --split-files "$accession" -O "$outdir"
    if [ $? -ne 0 ]; then
      echo "Ошибка при скачивании SRA: $accession"
      exit 1
    fi
    if [ -f "$outdir/${accession}_1.fastq" ]; then
        mv "$outdir/${accession}_1.fastq" "$outdir/${accession}_R1.fastq"
    fi
    if [ -f "$outdir/${accession}_2.fastq" ]; then
        mv "$outdir/${accession}_2.fastq" "$outdir/${accession}_R2.fastq"
    fi
}

#########################
# Разбор аргументов
input_files=()
outdir=""
concurrency=""

while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    -i)
      shift
      while [[ $# -gt 0 && "$1" != -* ]]; do
         input_files+=("$1")
         shift
      done
      ;;
    -o)
      outdir="$2"
      shift 2
      ;;
    -t)
      concurrency="$2"
      shift 2
      ;;
    *)
      echo "Неизвестный параметр: $1"
      exit 1
      ;;
  esac
done

if [[ -z "$outdir" || -z "$concurrency" || ${#input_files[@]} -eq 0 ]]; then
  echo "Usage: $0 -i <input_folder_or_files_or_txt> -o <output_directory> -t <concurrency>"
  exit 1
fi

#########################
# Определение режима работы и подготовка входных данных
# Режимы:
#   "sra"   – скачивание файлов по списку accession'ов;
#   "fastq" – работа с уже имеющимися FASTQ/FASTA файлами.
mode=""
if [[ ${#input_files[@]} -eq 1 && -f "${input_files[0]}" ]]; then
    first_line=$(head -n 1 "${input_files[0]}")
    # Если первая строка оканчивается на .fastq, .fq или .fasta (без учета регистра),
    # считаем, что передан список FASTQ/FASTA файлов.
    if [[ "$first_line" =~ \.(fastq|fq|fasta)$ ]]; then
         echo "Обработка текстового файла со списком FASTQ/FASTA файлов..."
         mode="fastq"
         mkdir -p "$outdir/FASTQ"
         while IFS= read -r file; do
             [ -z "$file" ] && continue
             cp "$file" "$outdir/FASTQ/"
         done < "${input_files[0]}"
         input_folder="$outdir/FASTQ"
    else
         echo "Обработка текстового файла со списком SRA-акцессий..."
         mode="sra"
         samples=()
         while IFS= read -r line; do
             [[ -z "$line" ]] && continue
             samples+=("$line")
         done < "${input_files[0]}"
         mkdir -p "$outdir/FASTQ"
         for sample in "${samples[@]}"; do
             if [ ! -f "$outdir/FASTQ/${sample}_R1.fastq" ] && [ ! -f "$outdir/FASTQ/${sample}.fastq" ]; then
                 download_sra "$sample" "$outdir/FASTQ" &
             fi
             # Ограничиваем число параллельных задач
             while (( $(jobs -rp | wc -l) >= concurrency )); do
                 sleep 1
             done
         done
         wait
         input_folder="$outdir/FASTQ"
    fi
elif [[ ${#input_files[@]} -eq 1 && -d "${input_files[0]}" ]]; then
    echo "Обработка входной папки с FASTQ/FASTA файлами..."
    mode="fastq"
    mkdir -p "$outdir/FASTQ"
    cp "${input_files[0]}"/*.[Ff]* "$outdir/FASTQ/"
    input_folder="$outdir/FASTQ"
elif [[ ${#input_files[@]} -gt 1 ]]; then
    echo "Обработка списка FASTQ/FASTA файлов..."
    mode="fastq"
    mkdir -p "$outdir/FASTQ"
    for f in "${input_files[@]}"; do
         cp "$f" "$outdir/FASTQ/"
    done
    input_folder="$outdir/FASTQ"
else
    echo "Неверный формат входных данных."
    exit 1
fi

# Приводим имена файлов в $outdir/FASTQ к стандартному виду (префикс SRR)
if [ "$mode" == "fastq" ]; then
  echo "Проверка и корректировка имён файлов в $outdir/FASTQ..."
  for file in "$outdir/FASTQ"/*.[Ff]*; do
      [ -e "$file" ] || continue
      base=$(basename "$file")
      if [[ ! $base =~ ^SRR ]]; then
           newbase="SRR${base}"
           echo "Переименование $base -> $newbase"
           mv "$outdir/FASTQ/$base" "$outdir/FASTQ/$newbase"
      fi
  done
fi

echo "Параметры:"
echo "  Режим: $mode"
echo "  Input folder: $input_folder"
echo "  Output: $outdir"
echo "  Concurrency: $concurrency"

#########################
# Создание необходимых директорий внутри outdir
mkdir -p "$outdir/FASTQ" "$outdir/ref" "$outdir/BAM" \
         "$outdir/stats/samtools_stats" "$outdir/stats/mosdepth" "$outdir/stats/fastqc" \
         "$outdir/stats/bcftools_stats" "$outdir/VCF" "$outdir/BCF"

#########################
# Скачивание и подготовка референса Mycobacterium tuberculosis
REFERENCE_ACC="NC_000962.3"
REFERENCE_FASTA="$outdir/ref/${REFERENCE_ACC}.fasta"
REF_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz"

if [ ! -f "$REFERENCE_FASTA" ]; then
  echo "Скачивание референса Mycobacterium tuberculosis ($REFERENCE_ACC)..."
  wget -q -O "$REFERENCE_FASTA.gz" "$REF_URL"
  if [ $? -ne 0 ]; then
    echo "Ошибка при скачивании референса через wget!"
    exit 1
  fi
  gunzip -f "$REFERENCE_FASTA.gz"
  if [ $? -ne 0 ]; then
    echo "Ошибка при распаковке референса!"
    exit 1
  fi
fi

if [ ! -f "${REFERENCE_FASTA}.fai" ]; then
  samtools faidx "$REFERENCE_FASTA"
fi

echo "Индексирование референса..."
bwa index "$REFERENCE_FASTA"
if [ $? -ne 0 ]; then
  echo "Ошибка при индексировании референса!"
  exit 1
fi

#########################
# Построение массива сэмплов
if [[ "$mode" == "sra" ]]; then
  echo "Найдено ${#samples[@]} сэмплов (SRA): ${samples[*]}"
else
  declare -a samples_fastq
  echo "Обработка входной папки с FASTQ/FASTA файлами..."
  echo "Файлы в $input_folder:"
  ls -1 "$input_folder"
  
  shopt -s nullglob
  for f in "$input_folder"/*.[Ff][Aa][Ss][Tt][Qq] "$input_folder"/*.[Ff][Qq] "$input_folder"/*.[Ff][Aa][Ss][Tt][Aa] ; do
    [ -e "$f" ] || continue
    base=$(basename "$f")
    sample=$(echo "$base" | sed -E 's/(_R[12])?(\.((fastq)|(fq)|(fasta)))$//')
    if [[ ! " ${samples_fastq[*]} " =~ " ${sample} " ]]; then
      samples_fastq+=("$sample")
    fi
  done
  shopt -u nullglob
  samples=("${samples_fastq[@]}")
  echo "Найдено ${#samples[@]} сэмплов (FASTQ/FASTA): ${samples[*]}"
fi

#########################
# Функция для вариантного анализа (создание VCF)
call_variants() {
  sample="$1"
  echo "Вариантный анализ для сэмпла $sample..."
  BAM_OUT="$outdir/BAM/${sample}.bam"
  BCF_DIR="$outdir/BCF"
  VCF_DIR="$outdir/VCF"
  PILEUP="$BCF_DIR/${sample}.pileup.bcf"
  CALLS="$BCF_DIR/${sample}.calls.bcf"
  VCF_OUT="$VCF_DIR/${sample}.vcf.gz"
  bcftools mpileup --threads "$CORES_BCFT" --min-MQ 30 --ignore-overlaps --max-depth 3000 \
    -f "$REFERENCE_FASTA" "$BAM_OUT" -Ou -o "$PILEUP"
  bcftools call --threads "$CORES_BCFT" --multiallelic-caller --ploidy 1 --variants-only \
    -Ou -o "$CALLS" "$PILEUP"
  bcftools view --threads "$CORES_BCFT" --include 'QUAL>20 && DP>10' -Oz -o "$VCF_OUT" "$CALLS"
  bcftools index "$VCF_OUT"
  bcftools stats --threads "$CORES_BCFT" "$VCF_OUT" > "$outdir/stats/bcftools_stats/${sample}.bcftools.txt"
}

#########################
# Функция для обработки одного сэмпла
process_sample() {
  sample="$1"
  echo "Обработка сэмпла: $sample"
  
  if [ -f "$outdir/FASTQ/${sample}_R1.fastq" ]; then
    R1="$outdir/FASTQ/${sample}_R1.fastq"
    R2="$outdir/FASTQ/${sample}_R2.fastq"
  elif [ -f "$outdir/FASTQ/${sample}_R1.fasta" ]; then
    R1="$outdir/FASTQ/${sample}_R1.fasta"
    R2="$outdir/FASTQ/${sample}_R2.fasta"
  elif [ -f "$outdir/FASTQ/${sample}.fastq" ]; then
    R1="$outdir/FASTQ/${sample}.fastq"
    R2=""
  elif [ -f "$outdir/FASTQ/${sample}.fasta" ]; then
    R1="$outdir/FASTQ/${sample}.fasta"
    R2=""
  elif [ -d "$input_folder" ]; then
    if [ -f "$input_folder/${sample}_R1.fastq" ]; then
      R1="$input_folder/${sample}_R1.fastq"
      R2="$input_folder/${sample}_R2.fastq"
    elif [ -f "$input_folder/${sample}_R1.fasta" ]; then
      R1="$input_folder/${sample}_R1.fasta"
      R2="$input_folder/${sample}_R2.fasta"
    elif [ -f "$input_folder/${sample}.fastq" ]; then
      R1="$input_folder/${sample}.fastq"
      R2=""
    elif [ -f "$input_folder/${sample}.fasta" ]; then
      R1="$input_folder/${sample}.fasta"
      R2=""
    fi
  fi

  if [ -z "$R1" ]; then
    echo "Входной файл для сэмпла $sample не найден. Пропуск."
    return
  fi

  echo "Запуск FastQC для сэмпла $sample (исходные чтения)..."
  if [ -n "$R2" ]; then
    fastqc -t "$CORES_FASTQC" "$R1" "$R2" -o "$outdir/stats/fastqc"
  else
    fastqc -t "$CORES_FASTQC" "$R1" -o "$outdir/stats/fastqc"
  fi

  echo "Выравнивание сэмпла $sample по референсу..."
  BAM_OUT="$outdir/BAM/${sample}.bam"
  if [ -n "$R2" ]; then
    bwa mem -t "$CORES_BWA" -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA" "$REFERENCE_FASTA" "$R1" "$R2" \
      | samtools view -bS - \
      | samtools sort -@ "$CORES_SAMTOOLS" -o "$BAM_OUT" -
  else
    bwa mem -t "$CORES_BWA" -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA" "$REFERENCE_FASTA" "$R1" \
      | samtools view -bS - \
      | samtools sort -@ "$CORES_SAMTOOLS" -o "$BAM_OUT" -
  fi

  if [ $? -ne 0 ]; then
    echo "Ошибка при выравнивании сэмпла $sample!"
    return
  fi

  samtools index "$BAM_OUT"

  echo "Запуск FastQC для BAM-файла сэмпла $sample..."
  fastqc -t "$CORES_FASTQC" "$BAM_OUT" -o "$outdir/stats/fastqc"

  echo "Запуск samtools stats для сэмпла $sample..."
  samtools stats "$BAM_OUT" > "$outdir/stats/samtools_stats/${sample}.samtools.txt"

  echo "Запуск mosdepth для сэмпла $sample..."
  mosdepth --threads "$CORES_MOSDEPTH" --fast-mode --no-per-base --use-median "$outdir/stats/mosdepth/${sample}" "$BAM_OUT"

  call_variants "$sample"
}

#########################
# Параллельная обработка сэмплов с ограничением по количеству одновременно работающих задач
for sample in "${samples[@]}"; do
  process_sample "$sample" &
  # Если число активных фоновых задач достигает лимита, ждем 1 секунду перед запуском следующей
  while (( $(jobs -rp | wc -l) >= concurrency )); do
      sleep 1
  done
done
wait

echo "Запуск MultiQC..."
multiqc "$outdir/stats" -o "$outdir/stats/multiqc"

echo "Пайплайн завершён."
