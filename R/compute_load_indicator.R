#' Computing the Pesticide Load Indicator for pesticide application data
#' The provided functions will compute the Pesticide Load Indicator (PLI) as described in Kudsk et al. (2018) for pesticide application data provided by the user.
#' Computing the PLI requires information on applied pesticides in a table format, as well as information on fate, ecotoxicity and human health properties of applied pesticide products, as provided in the Pesticide Properties Database (PPDB) of the University of Hertfordshire. See below for a detailed description.
#' The PLI can either be computed using user supplied information on pesticide properties or by automatically including the information based on the PPDB. Access to the PPDb requires a license - see http://sitem.herts.ac.uk/aeru/ppdb/.

#' @importFrom stats aggregate

required_columns_products <- c(
  "crop",
  "formula",
  "product",
  "reference.sum.risk.scores",
  "sum.risk.score"
)


prepare.products <- function(products)
{
    for (i in 4:5) {
        name = required_columns_products[i]
        # products[, name] <- as.numeric(unlist(products[,name]))
        products[,name] <- as.numeric(sub(", ", ".", unlist(products[,name]), fixed=TRUE))
    }
    products
}


optional_columns_products <- c("amount.applied", "standard.dosage")


#' @title Default load factors
#' @export
default.load.factors <- list(
    Load.Factor.SCI=20,
    Load.Factor.BCF=2.5,
    Load.Factor.SoilDT50=2.5,
    Load.Factor.Birds=1,
    Load.Factor.Mammals=1,
    Load.Factor.Fish=30,
    Load.Factor.Aquatic.Invertebrates=30,
    Load.Factor.Algae=3,
    Load.Factor.Aquatic.Plants=3,
    Load.Factor.Earthworms=2,
    Load.Factor.Bees=100,
    Load.Factor.Fish.Chronic=3,
    Load.Factor.Aquatic.Invertebrates.Chronic=3,
    Load.Factor.Earthworms.Chronic=2
)


required_columns_substances <- c(
  "substance",
  "product",
  "concentration",

  "SCI.Grow",
  "Reference.SCI.Grow",
  "Load.Factor.SCI",

  "BCF",
  "Reference.BCF",
  "Load.Factor.BCF",

  "SoilDT50",
  "Reference.SoilDT50",
  "Load.Factor.SoilDT50",

  "Birds.Acute.LD50.mg.kg",
  "Reference.Value.Birds",
  "Load.Factor.Birds",

  "Mammals.Acute.Oral.LD50.mg.kg.BW.day",
  "Reference.Value.Mammals",
  "Load.Factor.Mammals",

  "Fish.Acute.96hr.LC50.mg.l",
  "Reference.Value.Fish",
  "Load.Factor.Fish",

  "Aquatic.Invertebrates.Acute.48hr.EC50.mg.l",
  "Reference.Value.Aquatic.Invertebrates",
  "Load.Factor.Aquatic.Invertebrates",

  "Algae.Acute.72hr.EC50.Growth.mg.l",
  "Reference.Value.Algae",
  "Load.Factor.Algae",

  "Aquatic.Plants.Acute.7d.EC50.mg.l",
  "Reference.Value.Aquatic.Plants",
  "Load.Factor.Aquatic.Plants",

  "Earthworms.Acute.14d.LC50.mg.kg",
  "Reference.Value.Earthworms",
  "Load.Factor.Earthworms",

  "BeesLD50",
  "Reference.Value.Bees",
  "Load.Factor.Bees",

  "Fish.Chronic.21d.NOEC.mg.l.corrected",
  "Reference.Value.Fish.Chronic",
  "Load.Factor.Fish.Chronic",

  "water.phase.DT50.days",

  "Aquatic.Invertebrates.Chronic.21d.NOEC.mg.l.correted",
  "Reference.Value.Aquatic.Invertebrates.Chronic",
  "Load.Factor.Aquatic.Invertebrates.Chronic",

  "Earthworms.Chronic.14d.NOEC..Reproduction.mg.kg.corrected",
  "Reference.Value.Earthworms.Chronic",
  "Load.Factor.Earthworms.Chronic"

)



required_columns_substances_reduced <- c("CAS.number",
                                  "substance",
                                  "product",
                                  "concentration",
                                  "Load.Factor.SCI",
                                  "Load.Factor.BCF",
                                  "Load.Factor.SoilDT50",
                                  "Load.Factor.Birds",
                                  "Load.Factor.Mammals",
                                  "Load.Factor.Fish",
                                  "Load.Factor.Aquatic.Invertebrates",
                                  "Load.Factor.Algae",
                                  "Load.Factor.Aquatic.Plants",
                                  "Load.Factor.Earthworms",
                                  "Load.Factor.Bees",
                                  "Load.Factor.Fish.Chronic",
                                  "Load.Factor.Aquatic.Invertebrates.Chronic",
                                  "Load.Factor.Earthworms.Chronic")


required_columns_products_reduced <- c("product", "standard.dosage", "Year")



prepare.substances <- function(substances)
{
    for (i in 3:length(required_columns_substances)) {
        name = required_columns_substances[i]
        substances[,name] <- as.numeric(sub(", ", ".", unlist(substances[,name]), fixed=TRUE))
    }
    substances
}


#' @title Compute Pesticide Load Indicator with user supplied information on pesticide properties
#'
#' @param products Dataframe with raw pesticide application data.
#' @param substances Dataframe describing active ingredients of the applied pesticide products, including their ecotoxicity, fate and human health properties.
#' @return Dataframe with pesticide indicators for each pesticide application 
#' indicated in the products dataframe.
#' Computes Pesticide Load Indicator (L) and its subindicators:
#' The Human Health Load (HL), Ecotoxicity Load (TL) and Fate Load (FL).
#' If standard dosages are provided the Standard Treatment Index (STI) and
#' the Pesticide Load Index (LI=STI*L) are also computed.
#' @examples
#' \dontrun{
#' # load the dataframe containing the pesticide use data.
#' products_user <- products.load()
#' # load the (user-supplied) dataframe with detailed information on used pesticides.
#' substances_user <- substances.load()
#'
#' # Compute the Pesticide Load Indicator and its sub-indicators using the user supplied data.
#' indicators_user <- compute_pesticide_load_indicator(substances = substances_user,
#' products= products_user)
#' }
#' @export

compute_pesticide_load_indicator <- function(substances, products) {
  check_substance_column_names(substances)
  check_products_column_names(products)

  products <- prepare.products(products)
  substances <- prepare.substances(substances)

  substances <- compute_fate_load(substances)
  substances <- compute_toxity_load(substances)

  products <- compute_health_load(products)
  products <- compute_pesticide_load(products, substances)

  if (all(optional_columns_products %in% names(products))) {
    products <- compute_load_index(products)
  }

  return(products)
}


#' @title Compute Pesticide Load Indicator using the Pesticide Properties database
#'
#' @param products Dataframe with raw pesticide application data.
#' @param substances Dataframe describing active ingredients of the applied pesticide products and their CAS number.
#' @param folder Folder with exported xlsx files from PPDB containing information on active ingredient properties.
#' @return Dataframe with pesticide indicators for each pesticide application 
#' indicated in the products dataframe.
#' Computes Pesticide Load Indicator (L) and its subindicators:
#' The Human Health Load (HL), Ecotoxicity Load (TL) and Fate Load (FL).
#' @examples
#' \dontrun{
#' # load the dataframe containing the pesticide use data.
#' products_ppdb <- products.load()[,c("product","crop","standard.dosage","formula")]
#' # Add information on the year in which the product is used.
#' products_ppdb$Year <- c(2009,2010,2011,2012)
#'
#' # load the (user-supplied) dataframe with information on used pesticides
#' substances_ppdb <- substances.load()[,c("substance","product","concentration")]
#'
#' # Add the CAS number of each active ingredient to link to the Pesticide Properties database.
#' substances_ppdb$CAS.number <- c("94361-06-5","141517-21-7","111988-49-9","467-69-6",
#'                                 "1918-00-9","94-74-6","21087-64-9","142459-58-3")
#'
#' Load.factors <- c("Load.Factor.SCI","Load.Factor.BCF","Load.Factor.SoilDT50",
#'        "Load.Factor.Birds","Load.Factor.Mammals","Load.Factor.Fish",
#'        "Load.Factor.Aquatic.Invertebrates","Load.Factor.Algae","Load.Factor.Aquatic.Plants",
#'        "Load.Factor.Earthworms","Load.Factor.Bees","Load.Factor.Fish.Chronic",
#'        "Load.Factor.Aquatic.Invertebrates.Chronic","Load.Factor.Earthworms.Chronic")
#' # Add the Load factors as defined in the Danish Pesticide Load indicator.
#' # Alternatively supply own values for the Load factor.
#' for (i in 1:length(Load.factors)){
#'   substances_ppdb[,Load.factors[i]]<- rep(times=dim(substances_ppdb)[[1]],
#'                                           substances.load()[1,Load.factors[i]])
#' }
#' # Indicate the folder containing the "General","Fate", "Human" and "Ecotox" tables of the PPPDB.
#' # Excel files (under the exact same name, e.g. Human.xlsx) are required.
#' # Attention, a licensed access to the PPPDB (Lewis et al., 2016) is required.
#' # Note that the "Fate" table should include a column indicating the "SCI.GROW" values.
#' folder <- "path"
#'
#' # Compute the Pesticide Load Indicator and its sub-indicators using the user supplied data.
#' indicators_ppdb  <- compute_pesticide_load_indicator_ppdb(substances = substances_ppdb,
#' products= products_ppdb, folder=folder)
#' }
#' @export



compute_pesticide_load_indicator_ppdb <- function(substances, products, folder) {

  check_columns(substances, required_columns_substances_reduced, c(), "substances");
  check_columns(products, required_columns_products_reduced, c(), "products");

  tables <- match.ppdb(substances, products, folder)
  substances <- tables$substances
  products <- tables$products
  products <- products[, names(products) != "Year"]
  return(compute_pesticide_load_indicator(substances, products))
}

#' @title Compute Pesticide Load Indicator for a single application using the Pesticide Properties database
#'
#' @param folder Folder with exported xlsx files from PPDB containing information on active ingredient properties.
#' @param product Product name of the applied pesticide.
#' @param year Application year.
#' @param formula Load formulation factor for applied product.
#' @param substances List or vector of active ingredient names of the applied pesticide.
#' @param cas_numbers List or vector of CAS numbers of the respective active ingredients.
#' @param concentrations List or vector of product concentrations of the respective active ingredients.
#' @param ... overrides for default Load factors.
#' @return Dataframe with pesticide indicators for a single pesticide application
#' indicated by the user.
#' Computes Pesticide Load Indicator (L) and its subindicators:
#' The Human Health Load (HL), Ecotoxicity Load (TL) and Fate Load (FL).
#' @examples
#' \dontrun{
#' # Indicate the folder containing the "General","Fate", "Human" and "Ecotox" tables of the PPPDB.
#' # Excel files (under the exact same name, e.g. Human.xlsx) are required.
#' # Attention, a licensed access to the PPPDB (Lewis et al., 2016) is required.
#' # Note that the "Fate" table should include a column indicating the "SCI.GROW" values
#' folder <- "path"
#'
#' # Compute the Pesticide Load indicator and its subindicators for a single product.
#' # Allows to optionally alter the predefined Load.Values.
#' compute_pesticide_load_indicator_single(
#'   folder=folder,
#'   product="Agora",
#'   year=2009,
#'   formula=1.5,
#'   substances=c("Cyproconazol", "Trifloxystrobin"),
#'   cas_numbers=c("94361-06-5", "141517-21-7"),
#'   concentrations = c(0.08,0.1875)
#' )
#' }
#' @export

compute_pesticide_load_indicator_single <- function(folder, product, year, formula,
                                                    substances, cas_numbers,
                                                    concentrations, ...) {

    products <- data.frame(product=product,
                           Year=year,
                           crop="",
                           standard.dosage="",
                           formula=formula)

    substances <- data.frame(CAS.number=cas_numbers,
                             substance=substances,
                             concentration=concentrations,
                             product=product)
    overrides <- list(...)
    for (name in names(default.load.factors)) {
        value <- overrides[[name]]
        if (is.null(value))
            value <- default.load.factors[[name]]
        substances[[name]] <- value
    }
    compute_pesticide_load_indicator_ppdb(substances, products, folder)
}


#' @title Check if column names of substances dataframe are valid
#'
#' @description checks for valid colum names and stops execution if problems are detected
#' @param substances Dataframe describing active ingredients of the applied pesticide products.
#' @return No return value
#'
#' @export

check_substance_column_names <- function(substances)
{
    check_columns(substances,
                  required_columns_substances,
                  c(),
                  "substances")
}


#' @title Check if column names of applied pesticide products dataframe are valid
#'
#' @description checks for valid colum names and stops execution if problems are detected
#' @param products Dataframe with raw pesticide application data.
#' @return No return value
#'
#' @export

check_products_column_names <- function(products)
{
    check_columns(products,
                  required_columns_products,
                  optional_columns_products,
                  "products")
}

.is.superset <- function(a, b)
{
    length(intersect(a, b)) == length(intersect(b, b))
}



check_columns <- function(data_frame, required, optional, name) {
  found_columns <- names(data_frame)

  if (.is.superset(found_columns, required)) {
    return()
  }

  if (.is.superset(found_columns, union(required, optional))) {
    return()
  }

  missing <- setdiff(required, found_columns)
  missing_str <- paste(missing, collapse = ", ")

  if (length(missing_str) == 0) {
    missing_str <- "none"
  }

  message <- sprintf(
    "%s data frame is not valid. missing columns: %s", name, missing_str
  )

  stop(message)
}


compute_fate_load <- function(substances) {
  degradation <- (substances$SCI.Grow
    / substances$Reference.SCI.Grow
    * substances$Load.Factor.SCI)

  bioaccumulation <- (substances$BCF
    / substances$Reference.BCF
    * substances$Load.Factor.BCF)

  sci_growth_index <- (substances$SoilDT50
    / substances$Reference.SoilDT50
    * substances$Load.Factor.SoilDT50)

  substances$U <- degradation
  substances$B <- bioaccumulation
  substances$P <- sci_growth_index
  substances$Fate.Load.substances <- substances$U + substances$B + substances$P

  return(substances)
}


compute_toxity_load <- function(substances) {
  short_term_effect_birds <- (substances$Reference.Value.Birds
    / substances$Birds.Acute.LD50.mg.kg
    * substances$Load.Factor.Birds)

  substances$Fa <- ifelse(is.finite(short_term_effect_birds),
    short_term_effect_birds,
    0
  )

  short_term_effect_mammals <- (substances$Reference.Value.Mammals
    / substances$Mammals.Acute.Oral.LD50.mg.kg.BW.day
    * substances$Load.Factor.Mammals)

  substances$Pa <- ifelse(is.finite(short_term_effect_mammals),
    short_term_effect_mammals,
    0
  )

  short_term_effect_fish <- (substances$Reference.Value.Fish
    / substances$Fish.Acute.96hr.LC50.mg.l
    * substances$Load.Factor.Fish)

  substances$Fla <- ifelse(is.finite(short_term_effect_fish),
    short_term_effect_fish,
    0
  )

  short_term_effect_daphina <- (substances$Reference.Value.Aquatic.Invertebrates
    / substances$Aquatic.Invertebrates.Acute.48hr.EC50.mg.l
    * substances$Load.Factor.Aquatic.Invertebrates)

  substances$Da <- ifelse(is.finite(short_term_effect_daphina),
    short_term_effect_daphina,
    0
  )

  short_term_effect_algae <- (substances$Reference.Value.Algae
    / substances$Algae.Acute.72hr.EC50.Growth.mg.l
    * substances$Load.Factor.Algae)

  substances$Aa <- ifelse(is.finite(short_term_effect_algae),
    short_term_effect_algae,
    0
  )
  short_term_effect_aquatic_plants <- (substances$Reference.Value.Aquatic.Plants
    / substances$Aquatic.Plants.Acute.7d.EC50.mg.l
    * substances$Load.Factor.Aquatic.Plants)

  substances$Vp <- ifelse(is.finite(short_term_effect_aquatic_plants),
    short_term_effect_aquatic_plants,
    0
  )

  # short term effect Earthworms
  short_term_effect_earthworms <- (substances$Reference.Value.Earthworms
    / substances$Earthworms.Acute.14d.LC50.mg.kg
    * substances$Load.Factor.Earthworms)

  substances$Ra <- ifelse(is.finite(short_term_effect_earthworms),
    short_term_effect_earthworms,
    0
  )

  short_term_effect_bees <- (substances$Reference.Value.Bees
    / substances$BeesLD50
    * substances$Load.Factor.Bees)

  substances$Ba <- ifelse(is.finite(short_term_effect_bees),
    short_term_effect_bees,
    0
  )

  degradation_factor_water <- ((1 - exp((-log(2) / substances$water.phase.DT50.days) * 7))
                               / (((log(2) / substances$water.phase.DT50.days) * 7)))

  substances$Degradation.Factor.Water <- (
    ifelse(substances$water.phase.DT50.days == 0 | substances$water.phase.DT50.days == 708,
      1,
      degradation_factor_water
    )
  )

  long_term_effect_fish <- (substances$Reference.Value.Fish.Chronic
    / substances$Fish.Chronic.21d.NOEC.mg.l.corrected
    * substances$Load.Factor.Fish.Chronic
    * substances$Degradation.Factor.Water)

  substances$Flk <- ifelse(is.finite(long_term_effect_fish),
    long_term_effect_fish,
    0
  )

  long_term_effect_daphina <- (substances$Reference.Value.Aquatic.Invertebrates.Chronic
    / substances$Aquatic.Invertebrates.Chronic.21d.NOEC.mg.l.correted
    * substances$Load.Factor.Aquatic.Invertebrates.Chronic
    * substances$Degradation.Factor.Water)

  substances$Dk <- ifelse(is.finite(long_term_effect_daphina),
    long_term_effect_daphina,
    0
  )

  degradation_factor_soil <-
    (1 - exp((-log(2) / substances$SoilDT50) * 180)) / ((log(2) / substances$SoilDT50) * 180)

  substances$Degradation.Factor.Soil <- ifelse(substances$SoilDT50 == 0 | substances$SoilDT50 == 708,
    1,
    degradation_factor_soil
  )

  long_term_effect_earthworms <- (substances$Reference.Value.Earthworms.Chronic
    / substances$Earthworms.Chronic.14d.NOEC..Reproduction.mg.kg.corrected
    * substances$Load.Factor.Earthworms.Chronic
    * substances$Degradation.Factor.Soil)

  substances$Rk <- ifelse(is.finite(long_term_effect_earthworms),
    long_term_effect_earthworms,
    0
  )

  substances$Environmental.Toxicity.Substance <- (
    substances$Fa
      + substances$Pa
      + substances$Fla
      + substances$Da
      + substances$Aa
      + substances$Vp
      + substances$Ra
      + substances$Ba
      + substances$Flk
      + substances$Dk
      + substances$Rk
  )

  return(substances)
}


compute_health_load <- function(products) {
  products$HL <- (products$formula * products$sum.risk.score
    / products$reference.sum.risk.scores)

  return(products)
}


compute_pesticide_load <- function(products, substances) {

  TL.product <- substances$concentration * substances$Environmental.Toxicity.Substance
  FL.product <- substances$concentration * substances$Fate.Load.substances

  TL <- aggregate(TL.product, by = list(product = substances$product), FUN = sum)
  FL <- aggregate(FL.product, by = list(product = substances$product), FUN = sum)

  products$TL <- merge(products, TL, by="product", all.x=TRUE)$x
  products$FL <- merge(products, FL, by="product", all.x=TRUE)$x

  products$L <- products$HL + products$TL + products$FL

  return(products)
}


compute_load_index <- function(products) {
  sti_quotient <- products$amount.applied / products$standard.dosage
  load_index <- sti_quotient * products$L

  products$STI <- sti_quotient
  products$LI <- load_index

  return(products)
}

