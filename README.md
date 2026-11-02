# Домашнее задание 4. Клиническая аннотация вариантов

Вариант: `student_08.vcf` (50 SNV, GRCh38).

## Что лежит в репозитории

- `student_08.vcf` — исходный VCF.
- `results/annotated_variants_ru.tsv` — сводная таблица, 50 строк: идентификатор, координаты, ген, последствие VEP, правило выбора транскрипта, rsID, частота gnomAD, запись ClinVar и клиническая интерпретация.
- `results/annotated_variants.tsv` — та же таблица с HGVS, SIFT, PolyPhen, AlphaMissense, частотами экзомов и геномов и статусом ревью ClinVar.
- `results/vep_raw.json`, `results/vep_transcripts.tsv` — полный ответ Ensembl VEP, включая все транскрипты.
- `results/clinvar_matches.vcf`, `results/clinvar_annotation.tsv` — совпадения с ClinVar GRCh38.
- `results/gnomad_raw.json`, `results/gnomad_annotation.tsv` — gnomAD v4.
- `results/annotation_settings.md` — параметры запросов.
- `report.md` — отчёт: статистика VCF, разбор трёх вариантов и ограничения баз.
- `annotate.py` — скрипт, которым собраны таблицы. Повторный запуск использует уже сохранённые JSON VEP и gnomAD.

## Как повторить

Нужны Python 3, `tabix` и доступ к Ensembl REST, gnomAD и FTP NCBI.

```bash
python3 annotate.py
```

Отчёт написан по таблицам в `results/`. Перед отправкой в письмо преподавателю нужно подставить ФИО и ссылку на этот репозиторий.
# bioinfo4
# bioinfo4
