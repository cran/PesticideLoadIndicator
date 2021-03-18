#' @importFrom readxl read_excel
#' @importFrom stringr str_detect


read.excel <- function(excel_file) {
    table <- as.data.frame(read_excel(excel_file))
    names(table) <- make.names(names(table))
    table
}


extend.fate <- function(fate,ecotox)
{
    required_f <- c(
                  "LogP",
                  "SCI.GROW",
                  "Soil.DT50.typical...days",
                  "Soil.DT50.lab...days",
                  "Soil.DT50.notes",
                  "Water.phase.DT50...days",
                  "Active",
                  "ID"
                  )

    missing_f <- setdiff(required_f, names(fate))
    if (length(missing_f) > 0) {
        stop(paste("columns", paste(missing_f, collapse=", "), "missing in fate table"))
    }

     required_eco <- c(
                  "Bioconcentration.Factor",
                  "ID"
                  )

    missing_eco <- setdiff(required_eco, names(ecotox))
    if (length(missing_eco) > 0) {
        stop(paste("columns", paste(missing_eco, collapse=", "), "missing in ecotox table"))
    }

    fate$LogP[is.na(fate$LogP)] <- 0
    fate$LogP<-as.numeric(fate$LogP)
    fate$LogP[is.na(fate$LogP)] <- 0

    eco<- ecotox[,required_eco]
    fate<- merge(fate,eco,by= "ID")

    fate$BCF[fate$LogP >6 & fate$LogP !=0]<-10^((-0.2*fate$LogP[fate$LogP >6 & fate$LogP !=0])+(2.74*fate$LogP[fate$LogP >6 & fate$LogP !=0])-4.72)
    fate$BCF[fate$LogP <6 & fate$LogP !=0]<-10^((0.86*fate$LogP[fate$LogP <6 & fate$LogP !=0])-0.7)
    fate$BCF[fate$Bioconcentration.Factor!="" & !(is.na(fate$Bioconcentration.Factor))]<-fate$Bioconcentration.Factor[fate$Bioconcentration.Factor!="" & !(is.na(fate$Bioconcentration.Factor))]
    fate$BCF[is.na(fate$BCF)]<-0
    fate$BCF[fate$BCF=="Low risk"]<-0
    suppressWarnings(fate$BCF<-as.numeric(fate$BCF))
    fate$BCF[is.na(fate$BCF)]<-0
    fate$BCF[fate$BCF>5100]<-5100

    fate$SCI.GROW[fate$SCI.GROW=="Cannot be calculated"] <- 0
    fate$SCI.GROW[is.na(fate$SCI.GROW)] <- 0
    fate$SCI.GROW <- as.numeric(fate$SCI.GROW)
    fate$SCI.GROW[fate$SCI.GROW>1000] <- 0

    x<-cbind(as.numeric(fate$Soil.DT50.lab...days),as.numeric(fate$Soil.DT50.typical...days))
    fate$SoilDT50<-rowMeans(x,na.rm=TRUE)
    fate$stable<-ifelse(str_detect(fate$Soil.DT50.notes, "Stable"),1,0)
    fate$Stable<-ifelse(str_detect(fate$Soil.DT50.notes, "stable"),1,0)
    fate$stable<-fate$stable + fate$Stable
    fate$SoilDT50[fate$stable!="0"]<-2*354

    fate$SoilDT50[fate$SoilDT50==""] <- 0
    fate$SoilDT50[is.na(fate$SoilDT50)] <- 0
    fate$SoilDT50[fate$SoilDT50=="#N/A"] <- 0
    fate$SoilDT50[fate$SoilDT50>709] <- 0
    fate$SoilDT50 <- as.numeric(fate$SoilDT50)
    fate$SoilDT50[fate$Soil.DT50.notes=="Both iron and phosphate naturally occur in soil. Degradation will be very slow"] <- 0
    fate$SoilDT50[str_detect(fate$Active, "Copper")]<-0
    fate$SoilDT50[str_detect(fate$Active, "copper")]<-0
    fate$SoilDT50[str_detect(fate$Active, "Sulphur")]<-0
    fate$SoilDT50[str_detect(fate$Active, "sulphur")]<-0
    fate$SoilDT50[str_detect(fate$Active, "Iron")]<-0
    fate$SoilDT50[str_detect(fate$Active, "iron")]<-0



    fate$Water.phase.DT50...days[is.na(fate$Water.phase.DT50...days)]<-0
    fate$Water.phase.DT50...days[fate$Water.phase.DT50...days==""]<-0
    fate$Water.phase.DT50...days[fate$Water.phase.DT50...days=="<1"]<-0.5
    fate$Water.phase.DT50...days[fate$Water.phase.DT50...days==">100"]<-100
    fate$Water.phase.DT50...days[fate$Water.phase.DT50...days=="Slow, DT50 25-30 days"]<-30
    fate$Water.phase.DT50...days[fate$Water.phase.DT50...days=="Stable"]<-708

    fate

}



compute_R <- function(human) {
    # EC.Risk.Classification
    if (!('EC.Risk.Classification' %in% names(human))) {
        stop('Human table needs column EC.Risk.Classification')
    }
    R22 <- ifelse(str_detect(human$EC.Risk.Classification, "22"), 10, 0)
    R37 <- ifelse(str_detect(human$EC.Risk.Classification, "37"), 10, 0)
    R38 <- ifelse(str_detect(human$EC.Risk.Classification, "38"), 10, 0)
    R65 <- ifelse(str_detect(human$EC.Risk.Classification, "65"), 10, 0)
    R66 <- ifelse(str_detect(human$EC.Risk.Classification, "66"), 10, 0)
    R20 <- ifelse(str_detect(human$EC.Risk.Classification, "20"), 15, 0)
    R21 <- ifelse(str_detect(human$EC.Risk.Classification, "21"), 15, 0)
    R36 <- ifelse(str_detect(human$EC.Risk.Classification, "36"), 15, 0)
    R43 <- ifelse(str_detect(human$EC.Risk.Classification, "43"), 20, 0)
    R33 <- ifelse(str_detect(human$EC.Risk.Classification, "33"), 30, 0)
    R67 <- ifelse(str_detect(human$EC.Risk.Classification, "67"), 30, 0)
    R25 <- ifelse(str_detect(human$EC.Risk.Classification, "25"), 50, 0)
    R42 <- ifelse(str_detect(human$EC.Risk.Classification, "42"), 50, 0)
    R64 <- ifelse(str_detect(human$EC.Risk.Classification, "64"), 50, 0)
    R23 <- ifelse(str_detect(human$EC.Risk.Classification, "23"), 70, 0)
    R24 <- ifelse(str_detect(human$EC.Risk.Classification, "24"), 70, 0)
    R28 <- ifelse(str_detect(human$EC.Risk.Classification, "28"), 70, 0)
    R34 <- ifelse(str_detect(human$EC.Risk.Classification, "34"), 70, 0)
    R40 <- ifelse(str_detect(human$EC.Risk.Classification, "40"), 70, 0)

    R62 <- ifelse(str_detect(human$EC.Risk.Classification, "62"), 70, 0)
    R63 <- ifelse(str_detect(human$EC.Risk.Classification, "63"), 70, 0)
    R68 <- ifelse(str_detect(human$EC.Risk.Classification, "68"), 70, 0)
    R26 <- ifelse(str_detect(human$EC.Risk.Classification, "26"), 100, 0)
    R27 <- ifelse(str_detect(human$EC.Risk.Classification, "27"), 100, 0)
    R35 <- ifelse(str_detect(human$EC.Risk.Classification, "35"), 100, 0)
    R39 <- ifelse(str_detect(human$EC.Risk.Classification, "39"), 100, 0)
    R45 <- ifelse(str_detect(human$EC.Risk.Classification, "45"), 100, 0)
    R46 <- ifelse(str_detect(human$EC.Risk.Classification, "46"), 100, 0)
    R48 <- ifelse(str_detect(human$EC.Risk.Classification, "48"), 100, 0)
    R49 <- ifelse(str_detect(human$EC.Risk.Classification, "49"), 100, 0)
    R60 <- ifelse(str_detect(human$EC.Risk.Classification, "60"), 100, 0)
    R61 <- ifelse(str_detect(human$EC.Risk.Classification, "61"), 100, 0)

    R <- (R20 + R21 + R22 + R23 + R24 + R25 + R26 + R27 + R28 + R33 + R34
        + R35 + R36 + R37 + R38 + R39 + R40 + R42 + R43 + R45 + R46
        + R48 + R49 + R60 + R61 + R62 + R63 + R64 + R65 + R66 + R67 + R68)
    R
}

compute_HR <- function(human) {

    # CLP.classification.2013
    if (!('CLP.classification.2013' %in% names(human))) {
        stop('Human table needs column CLP.classification.2013')
    }

    #If there are two different risk points per H, the higher one has been taken
    #This happened for H300, H314, H330, H310)
    H302 <- ifelse(str_detect(human$CLP.classification.2013, "302"), 10, 0)
    H335 <- ifelse(str_detect(human$CLP.classification.2013, "335"), 10, 0)
    H315 <- ifelse(str_detect(human$CLP.classification.2013, "315"), 10, 0)
    H304 <- ifelse(str_detect(human$CLP.classification.2013, "304"), 10, 0)
    H066 <- ifelse(str_detect(human$CLP.classification.2013, "066"), 10, 0)
    H332 <- ifelse(str_detect(human$CLP.classification.2013, "332"), 15, 0)
    H312 <- ifelse(str_detect(human$CLP.classification.2013, "312"), 15, 0)
    H319 <- ifelse(str_detect(human$CLP.classification.2013, "319"), 15, 0)
    H317 <- ifelse(str_detect(human$CLP.classification.2013, "317"), 20, 0)
    H336 <- ifelse(str_detect(human$CLP.classification.2013, "336"), 30, 0)
    H301 <- ifelse(str_detect(human$CLP.classification.2013, "301"), 50, 0)
    H334 <- ifelse(str_detect(human$CLP.classification.2013, "334"), 50, 0)
    H362 <- ifelse(str_detect(human$CLP.classification.2013, "362"), 50, 0)
    H331 <- ifelse(str_detect(human$CLP.classification.2013, "331"), 70, 0)
    H311 <- ifelse(str_detect(human$CLP.classification.2013, "311"), 70, 0)
    H314 <- ifelse(str_detect(human$CLP.classification.2013, "314"), 100, 0)
    H351 <- ifelse(str_detect(human$CLP.classification.2013, "351"), 70, 0)
    H318 <- ifelse(str_detect(human$CLP.classification.2013, "318"), 70, 0)
    H373 <- ifelse(str_detect(human$CLP.classification.2013, "373"), 70, 0)
    H361 <- ifelse(str_detect(human$CLP.classification.2013, "361"), 70, 0)
    H371 <- ifelse(str_detect(human$CLP.classification.2013, "371"), 70, 0)
    H341 <- ifelse(str_detect(human$CLP.classification.2013, "341"), 70, 0)
    H330 <- ifelse(str_detect(human$CLP.classification.2013, "330"), 100, 0)
    H300 <- ifelse(str_detect(human$CLP.classification.2013, "300"), 85, 0)
    H310 <- ifelse(str_detect(human$CLP.classification.2013, "310"), 100, 0)
    H370 <- ifelse(str_detect(human$CLP.classification.2013, "370"), 100, 0)
    H350 <- ifelse(str_detect(human$CLP.classification.2013, "350"), 100, 0)
    H340 <- ifelse(str_detect(human$CLP.classification.2013, "340"), 100, 0)
    H372 <- ifelse(str_detect(human$CLP.classification.2013, "372"), 100, 0)
    H360 <- ifelse(str_detect(human$CLP.classification.2013, "360"), 100, 0)

    HR <- (H066 + H300 + H301 + H302 + H304 + H310 + H311 + H312 + H314 + H315
          + H317 + H318 + H319 + H330 + H331 + H332 + H334 + H335 + H336 + H340
          + H341 + H350 + H351 + H360 + H361 + H362 + H370 + H371 + H372 + H373)
    HR
}


extend.products.table <- function(products_table, substances_table, human, general)
{
    if (!('ID' %in% names(human))) {
        stop('Human table needs column ID')
    }

    if (!('CASS.RN' %in% names(general))) {
        stop('General table needs column "CASS RN"')
    }

    cas.index <- match("CASS.RN", names(general))

    missing.cas <- c()
    for (irow in 1:nrow(products_table)) {
        sum.risk.score <- 0.0
        products_row = products_table[irow,]
        substances_rows = substances_table[substances_table$product == products_row$product,]

        for (jrow in 1:nrow(substances_rows)) {
            substance_row = substances_rows[jrow,]
            CAS = substance_row$CAS.number
            substance = substance_row$substance
            match = general[which(general[,cas.index] == CAS),]

            if (nrow(match) == 0) {
                missing.cas <- c(missing.cas, CAS);
                if (length(missing.cas) < 11)
                    warning(paste("no entry for CAS", CAS, "\n"))
                if (length(missing.cas) == 11)
                    warning("supress missing CAS matches from now on\n\n")
                next
            }


            id <- match$ID
            human_row = human[which(human$ID == id),]
            if (nrow(human_row) == 0) {
                warning(paste("no entry for human risk for id", id, "\n"))
                next
            }
            if (products_row$Year <= 2012) {
                sum.risk.score <- c(sum.risk.score ,compute_R(human_row))
            }
            else {
                sum.risk.score <- c(sum.risk.score, compute_HR(human_row))
            }
        }
        products_table[irow, "sum.risk.score"] <- max(sum.risk.score, na.rm=T)
        products_table[irow, "reference.sum.risk.scores"] <- 350
    }

    # remove duplicates
    missing.cas <- union(missing.cas, missing.cas)
    if (length(missing) > 0) {
        txt <- paste(missing.cas, collapse=", ")
        warning(paste("\nthe CAS numbers", txt, "caused problems, please fix this\n\n"))
    }

    products_table
}


create.substances.table <- function(input_table, general, fate, ecotox) {

    if (!('CASS.RN' %in% names(general))) {
        stop('General table needs column "CASS RN"')
    }
    cas.index <- match("CASS.RN", names(general))

    required.fate <- c("ID", "SCI.GROW", "Water.phase.DT50...days","Active")

    missing <- setdiff(required.fate, names(fate))
    if (length(missing) > 0) {
        stop(paste("columns", paste(missing, collapse=", "), "missing in fate table"))
    }

    required.ecotox <- c("ID",
                         "Birds...Acute.LD50.mg.kg",
                         "Mammals...Acute.Oral.LD50.mg.kg.BW.day",
                         "Fish...Acute.96hr.LC50.mg.l",
                         "Aquatic.Invertebrates...Acute.48hr.EC50.mg.l",
                         "Algae...Acute.72hr.EC50.Growth.mg.l",
                         "Aquatic.Plants...Acute.7d.EC50.mg.l",
                         "Earthworms...Acute.14d.LC50.mg.kg",
                         "Honeybees...Contact.acute.48hr.LD50.ug.per.bee",
                         "Honeybees...Oral.Acute.48hr.LD50.ug.per.bee",
                         "Fish...Chronic.21d.NOEC.mg.l",
                         "Aquatic.Invertebrates...Chronic.21d.NOEC.mg.l",
                         "Earthworms...Chronic.NOEC..Reproduction.mg.kg",
                         "Bioconcentration.Factor")

    missing <- setdiff(required.ecotox, names(ecotox))
    if (length(missing) > 0) {
        stop(paste("columns", paste(missing, collapse=", "), "missing in ecotox table"))
    }

    names(input_table)<-make.names(names(input_table))

    required.input <- c("substance", "product", "concentration")

    missing <- setdiff(required.input, names(input_table))
    if (length(missing) > 0) {
        stop(paste("columns", paste(missing, collapse=", "), "missing in substances table"))
    }

    fate <- extend.fate(fate,ecotox)

    col_names <- required_columns_substances

    result <- data.frame(matrix(NA, ncol = length(col_names)),
                        stringsAsFactors=FALSE)
    names(result) <- col_names

    missing.cas <- c()
    missing.ecotox <- c()

    for (irow in 1:nrow(input_table)) {
        row = input_table[irow,]
        CAS = row$CAS.number
        substance = row$substance
        if (CAS == "" || substance == "")
            next

        match = general[which(general[, cas.index] == CAS),]

        if (nrow(match) == 0) {
            if (!(CAS %in% missing.cas)) {
                missing.cas <- c(missing.cas, CAS);
                if (length(missing.cas) < 11)
                    warning(paste("no entry for CAS", CAS, "\n"))
                if (length(missing.cas) == 11)
                    message("supress missing CAS matches from now on\n\n")
            }
            next
        }

        id <- match$ID

        fate_row <- fate[which(fate$ID == id),]
        ecotox_row <- ecotox[which(ecotox$ID == id),]
        for (name in required.ecotox)
            if (is.na(ecotox_row[name])) {
                if (!(CAS %in% missing.ecotox)) {
                    missing.ecotox <- c(missing.ecotox, CAS);
                    if (length(missing.ecotox) < 11)
                        warning(paste("entry", name, "for CAS", CAS, "is NAN in ecotox\n"))
                    if (length(missing.ecotox) == 11)
                        message("supress missing ecotox data from now on\n\n")
                }
                next
            }

        if (is.na(ecotox_row$Honeybees...Contact.acute.48hr.LD50.ug.per.bee)) {
            ecotox_honeybees <- ecotox_row$Honeybees...Oral.Acute.48hr.LD50.ug.per.bee
        }
        else
            ecotox_honeybees <-
            mean(c(as.numeric(ecotox_row$Honeybees...Contact.acute.48hr.LD50.ug.per.bee),
             as.numeric(ecotox_row$Honeybees...Oral.Acute.48hr.LD50.ug.per.bee)),na.rm=T)

        new_row <- c(
                row$substance,
                row$product,
                row$concentration,

                fate_row$SCI.GROW,
                12.5,
                20, # row$Load.Factor.SCI,

                fate_row$BCF,
                5100,
                2.5, # row$Load.Factor.BCF,

                fate_row$SoilDT50,
                354,
                2.5, # row$Load.Factor.SoilDT50,

                ecotox_row$Birds...Acute.LD50.mg.kg,
                10,
                1, # row$Load.Factor.Birds,

                ecotox_row$Mammals...Acute.Oral.LD50.mg.kg.BW.day,
                8,
                1, # row$Load.Factor.Mammals,

                ecotox_row$Fish...Acute.96hr.LC50.mg.l,
                0.00021,
                30, # row$Load.Factor.Fish,

                ecotox_row$Aquatic.Invertebrates...Acute.48hr.EC50.mg.l,
                0.00015,
                30, # row$Load.Factor.Aquatic.Invertebrates,

                ecotox_row$Algae...Acute.72hr.EC50.Growth.mg.l,
                0.00002,
                3, # row$Load.Factor.Algae,

                ecotox_row$Aquatic.Plants...Acute.7d.EC50.mg.l,
                0.00035,
                3, # row$Load.Factor.Aquatic.Plants,

                ecotox_row$Earthworms...Acute.14d.LC50.mg.kg,
                1.0,
                2, # row$Load.Factor.Earthworms,

                ecotox_honeybees,
                0.015,
                100, # row$Load.Factor.Bees,

                ecotox_row$Fish...Chronic.21d.NOEC.mg.l,
                0.00015,
                3, # row$Load.Factor.Fish.Chronic,

                fate_row$Water.phase.DT50...days,

                ecotox_row$Aquatic.Invertebrates...Chronic.21d.NOEC.mg.l,
                0.00015,
                3, # row$Load.Factor.Aquatic.Invertebrates.Chronic,

                ecotox_row$Earthworms...Chronic.NOEC..Reproduction.mg.kg,
                0.8,
                2 # row$Load.Factor.Earthworms.Chronic
                )

        result[nrow(result) + 1, ] <- new_row

    }
    problematic.cas <- union(missing.cas, missing.ecotox)
    if (length(problematic.cas) > 0) {
        txt <- paste(problematic.cas, collapse=", ")
        warning(paste("\nthe CAS numbers", txt, "caused problems, please fix this\n\n"))
    }

    result;
}


#' @title Expend tables with information on ecotoxicity, fate and human health properties from PPDB
#'
#' @param products Dataframe with raw pesticide application data.
#' @param substances Dataframe describing active ingredients of the applied pesticide products and their CAS number.
#' @param folder Folder with exported xlsx files from PPDB containing information on active ingredient properties.
#' @return names Lists with updated substance and product data frames.
#'
#' @export

match.ppdb <- function(substances, products, folder) {

    suppressWarnings({
        human <- read.excel(file.path(folder, "Human.xlsx"))
        general <- read.excel(file.path(folder, "General.xlsx"))
        fate <- read.excel(file.path(folder, "Fate.xlsx"))
        ecotox <- read.excel(file.path(folder, "Ecotox.xlsx"))
    })

    products <- extend.products.table(products, substances, human, general)
    substances <- create.substances.table(substances, general, fate, ecotox)

    return(list(products=products, substances=substances))
}
