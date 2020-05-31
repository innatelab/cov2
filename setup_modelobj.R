if (!exists("modelobj")) {
message("Automatically defining modelobj and quantobj")
#modelobj <- "protgroup"
#quantobj <- "protgroup"
modelobj <- "protregroup"
quantobj <- "pepmodstate"
}

modelobjs_df <- msdata[[str_c(modelobj, "s")]]
modelobj_idcol <- str_c(modelobj, "_id")
quantobj_idcol <- str_c(quantobj, "_id")
# FIXME should be done in prepare_data
modelobjs_df$object_id <- modelobjs_df[[modelobj_idcol]]
modelobjs_df$object_label <- modelobjs_df[[str_c(modelobj, "_label")]]

modelobj2protein.df <- msdata[[str_c("protein2", modelobj)]]
modelobj2protein.df$object_id <- modelobj2protein.df[[modelobj_idcol]]

modelobj_suffix <- dplyr::case_when(modelobj == "protgroup" ~ "_pg",
                             modelobj == "protregroup" ~ "_prg",
                             modelobj == "ptmgroup" ~ "_ptm",
                             TRUE ~ NA_character_)

global_labu_shift <- get(str_c("global_", quantobj, "_labu_shift"))
obj_labu_min <- get(str_c(quantobj, "_labu_min"))
instr_calib <- get(str_c("instr_calib_", quantobj))
