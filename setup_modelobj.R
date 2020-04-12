if (!exists("modelobj")) {
message("Automatically defining modelobj and quantobj")
#modelobj <- "protgroup"
#quantobj <- "protgroup"
modelobj <- "protregroup"
quantobj <- "pepmodstate"
}

modelobjs_df <- msdata[[str_c(modelobj, "s")]]
modelobj_idcol <- str_c(modelobj, "_id")
# FIXME should be done in prepare_data
modelobjs_df$object_id <- modelobjs_df[[modelobj_idcol]]
modelobjs_df$object_label <- modelobjs_df[[str_c(modelobj, "_label")]]

modelobj2protein.df <- msdata[[str_c("protein2", modelobj)]]
modelobj2protein.df$object_id <- modelobj2protein.df[[modelobj_idcol]]

modelobj_suffix <- case_when(modelobj == "protgroup" ~ "_pg",
                             modelobj == "protregroup" ~ "_prg",
                             TRUE ~ NA_character_)

global_labu_shift <- case_when(quantobj == "pepmodstate" ~ global_pepmodstate_labu_shift,
                               quantobj == "protgroup" ~ global_protgroup_labu_shift,
                               TRUE ~ NA_real_)
instr_calib <- if (quantobj == "pepmodstate") {
    instr_calib_pepmodstate
} else if (quantobj == "protgroup") {
    instr_calib_protgroup
} else {
    NULL
}
