require(dplyr) # mutate, filter, select, ...
require(stringr) # str_c
require(data.table) # fread, fwrite
require(tidyverse)
require(logging)

Medication_Mapping <- function(Medication_DF = Medications){
  loginfo("Creating medication mapping...")
  # Create Groupings
  Oral_antihistamines <- c("Brompheniramine", "Chlorpheniramine","Dexbrompheniramine", "Dexchlorpheniramine",
                           "Fexofenadine", "Terfenadine", "Carbinoxamine", "Clemastine",
                           "Dimenhydrinate", "Diphenhydramine", "Doxylamine", "Tripelennamine",
                           "Astemizole", "Desloratadine", "Loratadine", "Promethazine",
                           "Cetirizine", "Hydroxyzine", "Levocetirizine", "Azatadine", "Cyproheptadine")
  Intranasal_corticosteroids <- c("Beclomethasone diproprionate", "Budesonide", "Ciclesonide", "Dexamethasone",
                                  "Flunisolide", "Fluticasone", "Mometasone", "Triamcinolone")
  Intranasal_antihistamines <- c("Azelastine", "Azelastine/Fluticasone", "Olopatadine")
  Intranasal_decongestants <- c("Oxymetazoline", "Phenylephrine", "Tetrahydrozoline")
  Intranasal_anticholinergics <- c("Ipratropium")
  Intranasal_chromones <- c("Cromolyn")
  Oral_decongestants <- c("Phenylephrine", "Phenylpropanolamine", "Pseudoephedrine")
  Oral_antihistamine_decongestant <- c("Acrivastine/pseudoephedrine", "Azatadine/pseudoephedrine",
                                       "Brompheniramine/phenylpropanolamine", "Brompheniramine/pseudoephedrine",
                                       "Carbinoxamine/pseudoephedrine", "Cetirizine/pseudoephedrine",
                                       "Chlorpheniramine/phenindamine/phenylpropanolamine", "Chlorpheniramine/phenylephrine", 
                                       "Chlorpheniramine/phenylephrine/phenylpropanolamine/phenyltoloxamine", "Chlorpheniramine/phenylpropanolamine",
                                       "Chlorpheniramine/pseudoephedrine", "Clemastine/phenylpropanolamine",
                                       "Desloratadine/pseudoephedrine", "Dexbrompheniramine/pseudoephedrine",
                                       "Dexchlorpheniramine/pseudoephedrine", "Fexofenadine/pseudoephedrine",
                                       "Loratadine/pseudoephedrine", "Phenylephrine/promethazine",
                                       "Pseudoephedrine/terfenadine", "Pseudoephedrine/triprolidine")
  Mepolizumab <- c("Mepolizumab")
  Omalizumab <- c("Omalizumab")
  Oral_corticosteroids <- c("Prednisolone", "Prednisone")
  IV_steroids <- c("Methylprednisolone")
  Inhaled_corticosteroids <- c("Beclomethasone dipropionate", "Budesonide", "Ciclesonide", "Dexamethasone",
                            "Flunisolide", "Fluticasone", "Mometasone", "Triamcinolone")
  Inhaled_chromones <- c("Nedocromil", "Cromolyn")
  Inhaled_anticholinergics <- c("Umeclidinium", "Aclidinium", "Ipratropium", "Tiotropium")
  Antileukotrienes <- c("Montelukast", "Zafirlukast", "Zileuton")
  Long_acting_beta_agonists <- c("Arformoterol", "Formoterol", "Indacaterol", "Salmeterol")
  ICS_LABA <- c("Budesonide/formoterol", "Fluticasone/salmeterol", "Fluticasone/vilanterol", "Formoterol/mometasone")
  Epinephrine <- c("Epinephrine", "Epinephrine,racemic")
  Short_acting_beta_agonists <- c("Albuterol", "Bitolterol", "Isoetharine", "Isoproterenol",
                                  "Levalbuterol", "Metaproterenol", "Pirbuterol", "Terbutaline")
  SABA_anticholinergics <- c("Albuterol/ipratropium", "Ipratropium/albuterol")
  Xanthines <- c("Aminophylline", "Oxtriphylline", "Theophylline")
  ## Sort Medications
  Medication_DF <- Medication_DF %>% mutate(Medication_Group = NA, Medication_Name = NA)
  logdebug(str_c(Medication_DF %>% filter(!is.na(Medication_Name)) %>% nrow(), " out of ", nrow(Medication_DF), " filled"))
  # Oral_antihistamines
  loginfo("Processing Oral Antihistamines...")
  Medication_DF <- Medication_DF %>%
    mutate(Medication_Name = ifelse(is.na(Medication_Name) & grepl("^Bromphenir", Medication), "Brompheniramine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Chlor(( |-)trimeton|pheniramine( (2|4|E|s|t|m.*LMR)|-o|\\)))", Medication), "Chlorpheniramine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Dexbromphen", Medication), "Dexbrompheniramine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Dexchlorpheniramine", Medication), "Dexchlorpheniramine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Allegra.[^d]|Fenofexidine|^Fexofenadine((?![Pp]seudoephedrine).)*$|Mucinex Allergy", Medication, perl = TRUE), "Fexofenadine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Terfenadine|Seldane", Medication), "Terfenadine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Carbinoxamine", Medication), "Carbinoxamine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Clemastine|Tavist", Medication), "Clemastine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("^Dimenhydrinate|^Dramamine", Medication), "Dimenhydrinate", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Banophen|^Benadryl( |-)[^\\dU]|Diphenhist|^((Allergy|Sleep Aid|Unisom) \\(|Id-)*Diphenhydramine( |-|\\))((?!Acetaminophen|MGH|[Zz]inc|chemo).)*$|-Dryl", Medication, perl = TRUE), "Diphenhydramine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Doxylamine [Ss]uccinate|Unisom sleeptabs", Medication), "Doxylamine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Tripelennamine|Pyribenzamine", Medication), "Tripelennamine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Astemizole|Hismanal", Medication), "Astemizole", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Desloratadine|Clarinex( |-)[^d]", Medication), "Desloratadine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("^.*Loratadine((?![Pp]seudoephedrine).)*$|Alavert|Claritin( |-)[^dD]|Wal-ltin|Allerclear|Id-Alisertib", Medication, perl = TRUE), "Loratadine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("^.*Promethazine[^/]((?!MGH|Codeine|dm|Syr In).)*$|^.*Phenergan((?!codeine|MGH).)*$|Phenadoz|Promethegan", Medication, perl = TRUE), "Promethazine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("^.*Cetirizine((?![Pp]seudoephedrine).)*$|Zyrtec( |-)[^dDI]|Mycostatin|[Nn]ystatin(((?![Tt]opical).)*cream|.*(k unit|topical.*oncall))|Nystop.*Pad", Medication, perl = TRUE), "Cetirizine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Atarax|Hydroxyzine((?! up).)*$|Vistaril", Medication, perl = TRUE), "Hydroxyzine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("(l|L)evocetirizine", Medication), "Levocetirizine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Azatadine", Medication), "Azatadine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Cyproheptadine|Periactin", Medication), "Cyproheptadine", Medication_Name),
           Medication_Group = ifelse(is.na(Medication_Group) & (Medication_Name %in% Oral_antihistamines), "Oral Antihistamines", Medication_Group))
  logdebug(str_c(Medication_DF %>% filter(!is.na(Medication_Name)) %>% nrow(), " out of ", nrow(Medication_DF), " filled"))
  # Intranasal_corticosteroids
  loginfo("Processing Intranasal Corticosteroids...")
  Medication_DF <- Medication_DF %>%
    mutate(Medication_Name = ifelse(is.na(Medication_Name) & grepl("Beclomethasone ([Dd]ipropionate)*.*([Nn]asal|80.*[Aa]erosol|cream)|Qnasl|(Beco|Vance)nase((?!inhaler).)*$", Medication, perl = TRUE), "Beclomethasone diproprionate", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Budesonide.*[Nn]as(a)*l|Rhinocort", Medication), "Budesonide", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Ciclesonide.*[Nn]asal|Omnaris|Alvesco \\(c|Zetonna", Medication), "Ciclesonide", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("^Dexamethasone (12|20) mg/100ml", Medication), "Dexamethasone", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Flunisolide.*[Nn]asal|Nasalide", Medication), "Flunisolide", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Fluticasone (.*(Inhl|Nasl|Suspension)|propionate (0.05|16))|[Ff]lonase|Veramyst.*Suspension", Medication), "Fluticasone", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Mometasone.*([Nn]asal|furoate((?!%|[Oo]intment|[Pp]owder|[Cc]ream|lotion).)*$)|Nasonex", Medication, perl = TRUE), "Mometasone", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Triamcinolone ([Aa]cetonide.*([Nn]as(a)*l|762|tube)|diacetate|hexacetonide|nasal)|Nasacort.([Aa][Qq]|.*(Nasl|oncall|Aerosol)$)|Azmacort.*oncall$|Kenalog inj", Medication), "Triamcinolone", Medication_Name),
           Medication_Group = ifelse(is.na(Medication_Group) & (Medication_Name %in% Intranasal_corticosteroids), "Intranasal Corticosteroids", Medication_Group))
  logdebug(str_c(Medication_DF %>% filter(!is.na(Medication_Name)) %>% nrow(), " out of ", nrow(Medication_DF), " filled"))
  # Intranasal_antihistamines
  loginfo("Processing Intranasal Antihistamines...")
  Medication_DF <- Medication_DF %>%
    mutate(Medication_Name = ifelse(is.na(Medication_Name) & grepl("[Aa]stelin|Azelastine.(hcl|.*([Nn]as))", Medication), "Azelastine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Azelastine.*[Ff]luticasone|Dymista", Medication), "Azelastine/Fluticasone", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Olopatadine.*Nas", Medication), "Olopatadine", Medication_Name),
           Medication_Group = ifelse(is.na(Medication_Group) & (Medication_Name %in% Intranasal_antihistamines), "Intranasal Antihistamines", Medication_Group))
  logdebug(str_c(Medication_DF %>% filter(!is.na(Medication_Name)) %>% nrow(), " out of ", nrow(Medication_DF), " filled"))
  # Intranasal_decongestants
  loginfo("Processing Intranasal Decongestants...")
  Medication_DF <- Medication_DF %>%
    mutate(Medication_Name = ifelse(is.na(Medication_Name) & grepl("Oxymetazoline (hcl|0.05)|A(ne)*frin|Genasal", Medication), "Oxymetazoline", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("^Phenylephrine (.*(%.*([Nn]asal|[Ss]pray))|hcl (0|1%|2.5%.*LMR))|Neo-synephrine", Medication), "Phenylephrine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Tetrahydrozoline.*([Nn]asal|hcl)|Tyzine", Medication), "Tetrahydrozoline", Medication_Name),
           Medication_Group = ifelse(is.na(Medication_Group) & (Medication_Name %in% Intranasal_decongestants), "Intranasal Decongestants", Medication_Group))
  logdebug(str_c(Medication_DF %>% filter(!is.na(Medication_Name)) %>% nrow(), " out of ", nrow(Medication_DF), " filled"))
  # Intranasal_anticholinergics
  loginfo("Processing Intranasal Anticholinergics...")
  Medication_DF <- Medication_DF %>%
    mutate(Medication_Name = ifelse(is.na(Medication_Name) & grepl("Ipratropium .*([Nn]asal|[Ss]ray)|Atrovent ([Nn]as|Hfa Inhl)", Medication), "Ipratropium", Medication_Name),
           Medication_Group = ifelse(is.na(Medication_Group) & (Medication_Name %in% Intranasal_anticholinergics), "Intranasal Anticholinergics", Medication_Group))
  logdebug(str_c(Medication_DF %>% filter(!is.na(Medication_Name)) %>% nrow(), " out of ", nrow(Medication_DF), " filled"))
  # Intranasal_chromones
  loginfo("Processing Intranasal Chromones...")
  Medication_DF <- Medication_DF %>%
    mutate(Medication_Name = ifelse(is.na(Medication_Name) & grepl("Cromolyn.*[Nn]asal|Nasalcrom", Medication), "Cromolyn", Medication_Name),
           Medication_Group = ifelse(is.na(Medication_Group) & (Medication_Name %in% Intranasal_chromones), "Intranasal Chromones", Medication_Group))
  logdebug(str_c(Medication_DF %>% filter(!is.na(Medication_Name)) %>% nrow(), " out of ", nrow(Medication_DF), " filled"))
  # Oral_decongestants
  loginfo("Processing Oral Decongestants...")
  Medication_DF <- Medication_DF %>%
    mutate(Medication_Name = ifelse(is.na(Medication_Name) & grepl("Phenylephrine (- LMR|.*([Oo]ral|[Tt]ablet))", Medication), "Phenylephrine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Phenylpropanolamine", Medication), "Phenylpropanolamine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("^Pseudoephedrine[^/-]((?![Tt]riprolidine).)*$|Sudafed|^Nasal Decongestant", Medication, perl = TRUE), "Pseudoephedrine", Medication_Name),
           Medication_Group = ifelse(is.na(Medication_Group) & (Medication_Name %in% Oral_decongestants), "Oral Decongestants", Medication_Group))
  logdebug(str_c(Medication_DF %>% filter(!is.na(Medication_Name)) %>% nrow(), " out of ", nrow(Medication_DF), " filled"))
  # Oral_antihistamine_decongestant 
  loginfo("Processing Oral Antihistamines and Oral Decongestants...")
  Medication_DF <- Medication_DF %>%
    mutate(Medication_Name = ifelse(is.na(Medication_Name) & grepl("Acrivastine.*[Pp]seudoephedrine|Semprex", Medication), "Acrivastine/pseudoephedrine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Azatadine.*[Pp]seudoephedrine", Medication), "Azatadine/pseudoephedrine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Dimetapp", Medication), "Brompheniramine/phenylpropanolamine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Brompheniramine.*[Pp]seudoephedrine|Bromfed", Medication), "Brompheniramine/pseudoephedrine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Rondec", Medication), "Carbinoxamine/pseudoephedrine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Cetirizine.*[Pp]seudoephedrine|Zyrtec-D", Medication), "Cetirizine/pseudoephedrine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Nolamine", Medication), "Chlorpheniramine/phenindamine/phenylpropanolamine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Chlorpheniramine.*[Pp]henylephrine|Ru-tuss|Rynatan", Medication), "Chlorpheniramine/phenylephrine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Naldecon", Medication), "Chlorpheniramine/phenylephrine/phenylpropanolamine/phenyltoloxamine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Chlorpheniramine.*[Pp]henylprop|Ornade", Medication), "Chlorpheniramine/phenylpropanolamine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Chlorpheniramine.*[Pp]seudoephed|Cpm|Deconamine|Fedahist|Pepain|[Pp]seudoephedrine.*[Cc]hlorpheniramine", Medication), "Chlorpheniramine/pseudoephedrine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Clemastine.*[Pp]henylpropanolamine|Tavist.[Dd]", Medication), "Clemastine/phenylpropanolamine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("[Dd]esloratadine.*[Pp]seudoephedrine|Clarinex-[Dd]", Medication), "Desloratadine/pseudoephedrine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Drixoral", Medication), "Dexbrompheniramine/pseudoephedrine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Dexchlorpheniramine.*[Pp]seudoephedrine", Medication), "Dexchlorpheniramine/pseudoephedrine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("[Ff]exofenadine.*[Pp]seudoephedrine|[Aa]llegra.[Dd]", Medication), "Fexofenadine/pseudoephedrine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Loratadine.*[Pp]seudoephedrine|(Claritin|Wal-ltin).[Dd]", Medication), "Loratadine/pseudoephedrine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Phenylephrine.*[Pp]romethazine|[Pp]romethazine.*[Pp]henylephrine", Medication), "Phenylephrine/promethazine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Seldane.[Dd]", Medication), "Pseudoephedrine/terfenadine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Pseudoephedrine.*[Tt]riprolidine|[Tt]riprolidine.*[Pp]seudoephedrine|Actifed|Unifed", Medication), "Pseudoephedrine/triprolidine", Medication_Name),
           Medication_Group = ifelse(is.na(Medication_Group) & (Medication_Name %in% Oral_antihistamine_decongestant), "Oral Antihistamine and Oral Decongestant", Medication_Group))
  logdebug(str_c(Medication_DF %>% filter(!is.na(Medication_Name)) %>% nrow(), " out of ", nrow(Medication_DF), " filled"))
  # Mepolizumab
  loginfo("Processing Mepolizumab...")
  Medication_DF <- Medication_DF %>%
    mutate(Medication_Name = ifelse(is.na(Medication_Name) & grepl("mepolizumab", Medication), "Mepolizumab", Medication_Name),
           Medication_Group = ifelse(is.na(Medication_Group) & (Medication_Name %in% Mepolizumab), "Mepolizumab", Medication_Group))
  logdebug(str_c(Medication_DF %>% filter(!is.na(Medication_Name)) %>% nrow(), " out of ", nrow(Medication_DF), " filled"))
  # Omalizumab
  loginfo("Processing Omalizumab...")
  Medication_DF <- Medication_DF %>%
    mutate(Medication_Name = ifelse(is.na(Medication_Name) & grepl("[Oo]malizumab((?!MGH).)*$", Medication, perl = TRUE), "Omalizumab", Medication_Name),
           Medication_Group = ifelse(is.na(Medication_Group) & (Medication_Name %in% Omalizumab), "Omalizumab", Medication_Group))
  logdebug(str_c(Medication_DF %>% filter(!is.na(Medication_Name)) %>% nrow(), " out of ", nrow(Medication_DF), " filled"))
  # Oral_corticosteroids
  loginfo("Processing Oral Corticosteroids...")
  Medication_DF <- Medication_DF %>%
    mutate(Medication_Name = ifelse(is.na(Medication_Name) & grepl("Methotrexate.*BWH$|Orapred|Prednisolone((?![Oo][Pp][Hh][Tt]|[Ee]ye|MGH|per 5 mg|1%).)*$|Pred forte|Prelone", Medication, perl = TRUE), "Prednisolone", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Predni[Ss][Oo][Nn][Ee]((?!MGH|oOral|ir).)*$|Rayos", Medication, perl = TRUE), "Prednisone", Medication_Name),
           Medication_Group = ifelse(is.na(Medication_Group) & (Medication_Name %in% Oral_corticosteroids), "Oral Corticosteroids", Medication_Group))
  logdebug(str_c(Medication_DF %>% filter(!is.na(Medication_Name)) %>% nrow(), " out of ", nrow(Medication_DF), " filled"))
  # IV_steroids
  loginfo("Processing IV Steroids...")
  Medication_DF <- Medication_DF %>%
    mutate(Medication_Name = ifelse(is.na(Medication_Name) & grepl("^(Depo.|Solu){0,1}[Mm]edrol((?!\\(MGH\\)).)*$|^(Id-)*Methylpred((?!\\(MGH\\)|oral|oint|HC).)*$|Osm_Mix Mannitol", Medication, perl = TRUE), "Methylprednisolone", Medication_Name),
           Medication_Group = ifelse(is.na(Medication_Group) & (Medication_Name %in% IV_steroids), "IV Steroids", Medication_Group))
  logdebug(str_c(Medication_DF %>% filter(!is.na(Medication_Name)) %>% nrow(), " out of ", nrow(Medication_DF), " filled"))
  # Inhaled_corticosteroids
  loginfo("Processing Inhaled Corticosteroids...")
  Medication_DF <- Medication_DF %>%
    mutate(Medication_Name = ifelse(is.na(Medication_Name) & grepl("Beclomethasone([Dd]ipr)*((?![Nn]asal|80.*Aerosol).)*$|Beclovent|Beconase inhaler|Qvar|Vancenase inhaler|Vanceril", Medication, perl = TRUE), "Beclomethasone dipropionate", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Budesonide(/formoterol fumarate 80|((?![Nn]asal|respule|Rect|3mg|[Ff]ormoterol|Dr|inhaler-oncall|FDA|unit).)*$)|[Pp]ulmicort", Medication, perl = TRUE), "Budesonide", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Ciclesonide|Alvesco 160", Medication), "Ciclesonide", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("^Dexamethasone sod phospha($|te 4 mg)", Medication), "Dexamethasone", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Flunisolide((?![Nn]asal).)*$|Aerobid", Medication, perl = TRUE), "Flunisolide", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Arnuity|Flovent|^Fluticasone((?!%|[Ss]almeterol|[Vv]ilanterol|[Nn]as(a)*l|Inhl|2243|16gm s|diskus|Top).)*$|Id-Fluticasone|Veramyst((?![Nn]asal).)*$", Medication, perl = TRUE), "Fluticasone", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Asmanex|Mometasone.*(Powder|Hfa Aerosol)", Medication), "Mometasone", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Azmacort((?!oncall).)*$|Nasacort 55 mcg/inh|Triamcinolone( inhaler-oncall|.*adap)$", Medication, perl = TRUE), "Triamcinolone", Medication_Name),
           Medication_Group = ifelse(is.na(Medication_Group) & (Medication_Name %in% Inhaled_corticosteroids), "Inhaled Corticosteroids", Medication_Group))
  logdebug(str_c(Medication_DF %>% filter(!is.na(Medication_Name)) %>% nrow(), " out of ", nrow(Medication_DF), " filled"))
  # Inhaled_chromones
  loginfo("Processing Inhaled Chromones...")
  Medication_DF <- Medication_DF %>%
    mutate(Medication_Name = ifelse(is.na(Medication_Name) & grepl("Nedocromil sodium", Medication), "Nedrocromil", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Cromolyn((?!%|[Nn]asal|[Oo]ral).)*$|Intal", Medication, perl = TRUE), "Cromolyn", Medication_Name),
           Medication_Group = ifelse(is.na(Medication_Group) & (Medication_Name %in% Inhaled_chromones), "Inhaled Chromones", Medication_Group))
  logdebug(str_c(Medication_DF %>% filter(!is.na(Medication_Name)) %>% nrow(), " out of ", nrow(Medication_DF), " filled"))
  # Inhaled_anticholinergics
  loginfo("Processing Inhaled anticholinergics...")
  Medication_DF <- Medication_DF %>%
    mutate(Medication_Name = ifelse(is.na(Medication_Name) & grepl("Incruse|Umeclidinium((?!Vilanterol).)*$", Medication, perl = TRUE), "Umeclidinium", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Aclidinium|Tudorza", Medication), "Aclidinium", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("^Atrovent((?![Nn]as|Hfa Inhl|%).)*$|^[Id-]*Ipratropium[^-](albuterol.*neb |((?!albuterol|[Ss]pray|[Mm]ist|[Nn]asal|MGH|Dose).)*$)", Medication, perl = TRUE), "Ipratropium", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Spiriva|Tiotropium((?!Olodaterol).)*$", Medication, perl = TRUE), "Tiotropium", Medication_Name),
           Medication_Group = ifelse(is.na(Medication_Group) & (Medication_Name %in% Inhaled_anticholinergics), "Inhaled Anticholinergics", Medication_Group))
  logdebug(str_c(Medication_DF %>% filter(!is.na(Medication_Name)) %>% nrow(), " out of ", nrow(Medication_DF), " filled"))
  # Antileukotrienes
  loginfo("Processing Antileukotrienes...")
  Medication_DF <- Medication_DF %>%
    mutate(Medication_Name = ifelse(is.na(Medication_Name) & grepl("Montelukast|Singulair", Medication), "Montelukast", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Accolate|Zafirlukast", Medication), "Zafirlukast", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Zileuton|Zyflo", Medication), "Zileuton", Medication_Name),
           Medication_Group = ifelse(is.na(Medication_Group) & (Medication_Name %in% Antileukotrienes), "Antileukotrienes", Medication_Group))
  logdebug(str_c(Medication_DF %>% filter(!is.na(Medication_Name)) %>% nrow(), " out of ", nrow(Medication_DF), " filled"))
  # Long_acting_beta_agonists
  loginfo("Processing Long Acting Beta Agonists...")
  Medication_DF <- Medication_DF %>%
    mutate(Medication_Name = ifelse(is.na(Medication_Name) & grepl("Arformoterol|Brovana", Medication), "Arformoterol", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("^Formoterol|Foradil", Medication), "Formoterol", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Indacaterol", Medication), "Indacaterol", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("^Salmet|Serevent((?!21).)*$", Medication, perl = TRUE), "Salmeterol", Medication_Name),
           Medication_Group = ifelse(is.na(Medication_Group) & (Medication_Name %in% Long_acting_beta_agonists), "Long Acting Beta Agonists", Medication_Group))
  logdebug(str_c(Medication_DF %>% filter(!is.na(Medication_Name)) %>% nrow(), " out of ", nrow(Medication_DF), " filled"))
  # ICS_LABA
  loginfo("Processing Inhaled Corticosteroids and Long Acting Beta Agonists...")
  Medication_DF <- Medication_DF %>%
    mutate(Medication_Name = ifelse(is.na(Medication_Name) & grepl("[Bb]udesonide.*[Ff]ormoterol((?!6.9).)*$|Symbicort", Medication, perl = TRUE), "Budesonide/formoterol", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Advair|[Ff]luticasone.*[Ss]almeterol", Medication), "Fluticasone/salmeterol", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Breo|[Ff]luticasone.*[Vv]ilanterol", Medication), "Fluticasone/vilanterol", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Dulera|[Mm]ometasone.*[Ff]ormoterol", Medication), "Formoterol/mometasone", Medication_Name),
           Medication_Group = ifelse(is.na(Medication_Group) & (Medication_Name %in% ICS_LABA), "Inhaled Corticosteroids and Long Acting Beta Agonists", Medication_Group))
  logdebug(str_c(Medication_DF %>% filter(!is.na(Medication_Name)) %>% nrow(), " out of ", nrow(Medication_DF), " filled"))
  # Epinephrine
  loginfo("Processing Epinephrine...")
  Medication_DF <- Medication_DF %>%
    mutate(Medication_Name = ifelse(is.na(Medication_Name) & grepl("^Epinephrine.(1(.*inhalation| to 1)|B|oncall$)", Medication), "Epinephrine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Racepinephrin", Medication), "Epinephrine,racemic", Medication_Name),
           Medication_Group = ifelse(is.na(Medication_Group) & (Medication_Name %in% Epinephrine), "Epinephrine", Medication_Group))
  logdebug(str_c(Medication_DF %>% filter(!is.na(Medication_Name)) %>% nrow(), " out of ", nrow(Medication_DF), " filled"))
  # Short_acting_beta_agonists
  loginfo("Processing Short Acting Beta Agonists...")
  Medication_DF <- Medication_DF %>%
    mutate(Medication_Name = ifelse(is.na(Medication_Name) & grepl("^(Id-)*Albuterol[^,]((?!MGH|LMR (25|3846)|[Ii]pra|neb$|17gm can|extend|oncall|non-com|%ML).)*$|Proair|Proventil|Ventolin", Medication, perl = TRUE), "Albuterol", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Bitolterol", Medication), "Bitolterol", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Bronkometer|Isoetharine", Medication), "Isoetharine", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Isoproterenol hcl (-|1m)", Medication), "Isoproterenol", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("^Levalbuterol((?![Cc]o).)*$|Xopenex", Medication, perl = TRUE), "Levalbuterol", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Alupent|Metaproterenol", Medication), "Metaproterenol", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Maxair|Pirbuterol", Medication), "Pirbuterol", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Brethine|^Terbutaline((?!Sub|MGH|sulfate 1 *mg).)*$", Medication, perl = TRUE), "Terbutaline", Medication_Name),
           Medication_Group = ifelse(is.na(Medication_Group) & (Medication_Name %in% Short_acting_beta_agonists), "Short Acting Beta Agonists", Medication_Group))
  logdebug(str_c(Medication_DF %>% filter(!is.na(Medication_Name)) %>% nrow(), " out of ", nrow(Medication_DF), " filled"))
  # SABA_anticholinergics
  loginfo("Processing Short Acting Beta Agonists and Anticholinergics...")
  Medication_DF <- Medication_DF %>%
    mutate(Medication_Name = ifelse(is.na(Medication_Name) & grepl("Albuterol((?!sulfate).)*[Ii]pratropium|Combivent|Duoneb|^Ipratropium.*[Aa]lbuterol((?!ml$|vial).)*$", Medication, perl = TRUE), "Albuterol/ipratropium", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Ipratropium/albuterol 3 ml", Medication), "Ipratropium/albuterol", Medication_Name),
           Medication_Group = ifelse(is.na(Medication_Group) & (Medication_Name %in% SABA_anticholinergics), "Short Acting Beta Agonists and Anticholinergics", Medication_Group))
  logdebug(str_c(Medication_DF %>% filter(!is.na(Medication_Name)) %>% nrow(), " out of ", nrow(Medication_DF), " filled"))
  # Xanthines
  loginfo("Processing Xanthines...")
  Medication_DF <- Medication_DF %>%
    mutate(Medication_Name = ifelse(is.na(Medication_Name) & grepl("Aminophylline", Medication), "Aminophylline", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Choledyl|Oxtriphylline", Medication), "Oxtriphylline", Medication_Name),
           Medication_Name = ifelse(is.na(Medication_Name) & grepl("Slo |Theo|Uniphyl", Medication), "Theophylline", Medication_Name),
           Medication_Group = ifelse(is.na(Medication_Group) & (Medication_Name %in% Xanthines), "Xanthines", Medication_Group))
  logdebug(str_c(Medication_DF %>% filter(!is.na(Medication_Name)) %>% nrow(), " out of ", nrow(Medication_DF), " filled"))
  loginfo("Mapping complete")
  return(Medication_DF)
}

Create_ACTH_Cortisol_DHEA_Output_Columns <- function(ACTH_Cortisol_DHEA_Group, Group_Header, DF_to_fill){
  # Splitting up Output Columns step into substeps because some median date selection doesn't work well with summarise
  OC_pregroup <- ACTH_Cortisol_DHEA_Group %>% group_by(EMPI, Reference_Units, Seq_Date, Seq_Time) %>%
    summarise(nResults_Per_Time = n(),
              MinResult = min(Result),
              MaxResult = max(Result),
              MedianResult = median(Result),
              MeanResult = mean(Result)) %>%
    arrange(Seq_Time) %>% group_by(EMPI, Reference_Units, Seq_Date) %>%
    summarise(nResults_Per_Date = sum(nResults_Per_Time),
              nTimes_Per_Date = n(),
              FirstMin = first(MinResult),
              FirstMax = first(MaxResult),
              FirstMedian = first(MedianResult),
              FirstMean = first(MeanResult),
              FirstTime = first(Seq_Time)) %>%
    mutate(Seq_Date_Time = ifelse(is.na(FirstTime), as.character(Seq_Date), str_c(Seq_Date, " ", FirstTime))) %>% 
    group_by(EMPI, Reference_Units) %>%
    rename(!!as.symbol(str_c(Group_Header, "_Reference_Units")) := Reference_Units)
  DF_to_fill <- left_join(DF_to_fill, OC_pregroup %>% summarise(), by = "EMPI")
  OC_pregroup <- OC_pregroup %>% group_by(EMPI) %>% arrange(Seq_Date_Time)
  
  Output_Columns <- OC_pregroup %>%
    summarise(!!as.symbol(str_c(Group_Header, "_nTotalDates")) := n(),
              !!as.symbol(str_c(Group_Header, "_nTotalDatesTimes")) := sum(nTimes_Per_Date),
              !!as.symbol(str_c(Group_Header, "_nTotalResults")) := sum(nResults_Per_Date),
              !!as.symbol(str_c(Group_Header, "_All_Seq_Date_Times")) := paste(Seq_Date_Time, collapse = ";"))
  DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
  Output_Columns <- OC_pregroup %>%
    summarise(!!as.symbol(str_c(Group_Header, "_Overall_Min_Result")) := min(FirstMin),
              !!as.symbol(str_c(Group_Header, "_Overall_Min_Result_Date_First")) := Seq_Date_Time[first(which(FirstMin == min(FirstMin)))],
              !!as.symbol(str_c(Group_Header, "_Overall_Min_Result_Date_Last")) := Seq_Date_Time[last(which(FirstMin == min(FirstMin)))],
              !!as.symbol(str_c(Group_Header, "_All_Min_Results")) := paste(FirstMin, collapse = ";"))
  DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
  Output_Columns <- OC_pregroup %>%
    summarise(!!as.symbol(str_c(Group_Header, "_Overall_Max_Result")) := max(FirstMax),
              !!as.symbol(str_c(Group_Header, "_Overall_Max_Result_Date_First")) := Seq_Date_Time[first(which(FirstMax == max(FirstMax)))],
              !!as.symbol(str_c(Group_Header, "_Overall_Max_Result_Date_Last")) := Seq_Date_Time[last(which(FirstMax == max(FirstMax)))],
              !!as.symbol(str_c(Group_Header, "_All_Max_Results")) := paste(FirstMax, collapse = ";"))
  DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
  Output_Columns <- OC_pregroup %>% filter(abs(median(FirstMedian) - FirstMedian) == median(FirstMedian) - FirstMedian) %>%
    arrange(EMPI, FirstMedian) %>%
    summarise(!!as.symbol(str_c(Group_Header, "_Overall_Median_Result")) := max(FirstMedian),
              !!as.symbol(str_c(Group_Header, "_Overall_Median_Result_Date_First_or_closest_below")) := Seq_Date_Time[first(which(FirstMedian == max(FirstMedian)))],
              !!as.symbol(str_c(Group_Header, "_Overall_Median_Result_Date_Last_or_closest_below")) := Seq_Date_Time[last(which(FirstMedian == max(FirstMedian)))])
  DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
  # All Median is done separately because the previous computation dumps rows that we don't want dumped for the paste
  Output_Columns <- OC_pregroup %>%
    summarise(!!as.symbol(str_c(Group_Header, "_All_Median_Results")) := paste(FirstMedian, collapse = ";"),
              !!as.symbol(str_c(Group_Header, "_Overall_Mean_Result")) := mean(FirstMean),
              !!as.symbol(str_c(Group_Header, "_All_Mean_Results")) := paste(FirstMean, collapse = ";"))
  
  DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
  loginfo(str_c(nrow(Output_Columns), " subjects have a ", Group_Header, " test performed with useful information"))
  rm(OC_pregroup, Output_Columns)
  return(DF_to_fill)
}