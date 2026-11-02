# Настройки аннотации

Дата запросов: 2026-11-02. Вход: `student_08.vcf`, GRCh38, 50 SNV. Генотипа образца в файле нет.

## Ensembl VEP

- Эндпоинт: `POST https://rest.ensembl.org/vep/homo_sapiens/region`
- Релиз: 116 (`results/ensembl_software.json`)
- Входная строка на вариант: `CHROM POS ID REF ALT . . .`
- Параметры: `canonical=1`, `mane=1`, `hgvs=1`, `numbers=1`, `protein=1`, `appris=1`, `ccds=1`, `uniprot=1`, `domains=1`, `xref_refseq=1`, `variant_class=1`, `AlphaMissense=1`, `flag_pick=1`
- `pick_order`: `mane_select,mane_plus_clinical,canonical,appris,tsl,biotype,ccds,rank,length`
- SIFT и PolyPhen возвращаются сервисом для миссенс-вариантов без отдельного флага
- Правило строки сводной таблицы: один транскрипт по флагу PICK; на нём сообщается самый тяжёлый термин Sequence Ontology. Самое тяжёлое последствие по всем транскриптам хранится отдельно
- Сырой ответ: `results/vep_raw.json`
- Все транскрипты: `results/vep_transcripts.tsv`

## ClinVar

- Файл: `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz`
- Дата файла в заголовке: 2026-10-25
- Доступ: `tabix` по интервалу `chrom:pos-pos`, затем точное совпадение REF и ALT
- Интерпретация: `CLNSIG` (агрегатная герминальная классификация). Конфликт отправителей: `CLNSIGCONF`
- Если позиция есть в VCF, а `CLNSIG` пуст, в таблице стоит «агрегатная классификация отсутствует»
- Совпавшие строки с заголовком: `results/clinvar_matches.vcf`
- Таблица: `results/clinvar_annotation.tsv`

## gnomAD

- Эндпоинт: `https://gnomad.broadinstitute.org/api`
- Набор: `gnomad_r4`, GRCh38
- Идентификатор варианта: `chrom-pos-ref-alt`
- Частота в сводной таблице: joint AC/AN (экзомы и геномы вместе)
- Отдельно сохранены частоты экзомов и геномов, фильтры QC и FAF95 popmax
- Фильтр `AC0` означает, что аллель не оставлен в отфильтрованном наборе; для такого варианта частота в сводную таблицу не подставляется
- Сырой ответ: `results/gnomad_raw.json`
- Таблица: `results/gnomad_annotation.tsv`

## Сводная таблица

- `results/annotated_variants.tsv` и `.csv` — 50 строк, все поля
- `results/annotated_variants_ru.tsv` — те же 50 вариантов, поля из условия задания
- Статистика VCF: `results/vcf_statistics.json`
- Машиночитаемая копия этих настроек: `results/annotation_settings.json`
