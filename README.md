---
title: "Copernicus PPI (10 m) – Double-Logistic Processing"
author: "[Philippe Choler]"
date: "Dernière mise à jour : 15/05/2025"
output:
  html_document:
    toc: true
    toc_float: true
    number_sections: true
    theme: readable
    css: styles.css
  pdf_document:
    toc: true
    number_sections: true
    latex_engine: xelatex
    extra_dependencies: ["georgia"]
fontsize: 12pt
mainfont: Georgia
geometry: margin=1in
---

<style>
body {font-family: Georgia, serif; text-align: justify;}
h1,h2,h3,h4,h5,h6 {font-family: Georgia, serif;}
</style>

## Présentation

Ce dépôt fournit un workflow **R** complet pour :

* Télécharger les **produits Copernicus HR-VPP 10 m** (Plant-Phenology Index) sur n’importe quelle emprise et plage d’années.  
* Assembler & reprojeter les **10 paramètres PPI**, puis estimer les 4 paramètres manquants du modèle logistique double (*ONSET, GROWTH, OFFSET, SENESC*).  
* Calculer :
  * la pile quotidienne de l’**IRG** (*Instantaneous Rate of Green-up*) ;
  * les **phases phénologiques** (1 Croissance, 2 Plateau, 3 Dépérissement, 4 Sénescence) ;
  * le raster **IRG-max** (valeur maximale et date).

Les sorties sont enregistrées dans `output/data_<site>/IRG/`.

---

## Bibliographie (sélection)

| Année | Référence | Contribution |
|-------|-----------|--------------|
| 2003 | Zhang *et al.* — RSE 84 | DL 6 paramètres (fondement HR-VPP) |
| 2006 | Beck *et al.* — RSE 100 | DL vs Fourier / Gauss |
| 2007 | Fisher & Mustard — RSE 109 | Validation MODIS-DL |
| 2012 | Elmore *et al.* — GCB 18 | Variante DL 7 paramètres |
| 2018 | Jönsson *et al.* — Remote Sens. 10 | Implémentation VITO HR-VPP |

---

## Arborescence du projet

```
📂 Copernicus_PPI_DoubleLogistic/
│-- 📜 config.R
│-- 📂 functions/
│   │-- Function_processing.R
│   │-- Function_estimate_parameters.R
│   └-- Function_calcul_indicators.R
│-- 📂 input/
│   └-- 📂 data_<site>/
│      │-- 📂 downloads/        # tuiles HR-VPP brutes (LAEA)
│      │-- 📂 extent/
│         └-- Masque_<site>.shp
└-- 📂 output/
   └-- 📂 data_<site>/
      │-- 📂 raw_data/          # PPI param reprojetés
      └-- 📂 IRG/
         ├-- IRG_season_<site>_<year>.tif
         ├-- PhenologyPhase_<site>_<year>.tif
         └-- IRG_max_<site>_<year>.tif
```
