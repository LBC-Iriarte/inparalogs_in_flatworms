library(edgeR)

setwd("C:/Users/El3ti/Documents/Trabajo_Ciencia/Mauri_Imparalogos/Schistosoma")

data <- gsub("_all_counts.tab", "" ,list.files(path=".", pattern="*_all_counts.tab"))

counts_file <- gsub(" ","", paste(data ,"_all_counts.tab"))
group_file <- gsub(" ","", paste(data ,"_grupos.tab"))

print_correlation_tab <- gsub(" ","", paste(data ,"_dif_expression.tab"))

# 1) Cargar los datos a R:
raw_counts <- read.table(counts_file, header = TRUE, sep = "\t", row.names="Term_Count")
group_list <- read.table(group_file, header = TRUE, sep = "\t")

#Groups <- c(group_list$Grupo)
#Groups

Groups <- c("A", "A", "A", "A", "A", "A", "A", "B", "B", "B", "C", "C", "C", "D", "D", "D", "E", "E", "E", "F", "F", "F")

# 2) Producción del objeto a utilizar por edgeR
edgeR_Input <- DGEList(counts=raw_counts, group=Groups)

# 3) Filtrado de la Librería:
# filterByExpr.default <- function(y, design=NULL, group=NULL, lib.size=NULL, min.count=10, min.total.count=15, ...)
keep<- filterByExpr(edgeR_Input, group=Groups)
edgeR_Input <- edgeR_Input[keep, , keep.lib.sizes=FALSE]

# 4) Normalizado
edgeR_Input<-calcNormFactors(edgeR_Input)

# 5) Diseño de matriz experimental
design_matrix <- model.matrix(~0+Groups, data=edgeR_Input$samples) # No termino de entender como funciona...
colnames(design_matrix) <- levels(edgeR_Input$samples$group)

# 6) Calculos de expresión diferencial
edgeR_Input <- estimateDisp(edgeR_Input,design_matrix)

# 7) Ajuste de la distribución binomial
fit <- glmQLFit(edgeR_Input,design_matrix)

# 8) Expresión diferencial:
# A	B
current_contrast <- makeContrasts(A-B, levels=design_matrix)
current_Test <- glmQLFTest(fit, contrast=current_contrast)
current_results <- topTags(current_Test,20000) 
write.table(current_results, file="Results_A_VS_B_2023-02-23.tab", sep = "\t", col.names = NA)

# A	C
current_contrast <- makeContrasts(A-C, levels=design_matrix)
current_Test <- glmQLFTest(fit, contrast=current_contrast)
current_results <- topTags(current_Test,20000) 
write.table(current_results, file="Results_A_VS_C_2023-02-23.tab", sep = "\t", col.names = NA)

# A	D
current_contrast <- makeContrasts(A-D, levels=design_matrix)
current_Test <- glmQLFTest(fit, contrast=current_contrast)
current_results <- topTags(current_Test,20000) 
write.table(current_results, file="Results_A_VS_D_2023-02-23.tab", sep = "\t", col.names = NA)

# A	E
current_contrast <- makeContrasts(A-E, levels=design_matrix)
current_Test <- glmQLFTest(fit, contrast=current_contrast)
current_results <- topTags(current_Test,20000) 
write.table(current_results, file="Results_A_VS_E_2023-02-23.tab", sep = "\t", col.names = NA)

# A	F
current_contrast <- makeContrasts(A-F, levels=design_matrix)
current_Test <- glmQLFTest(fit, contrast=current_contrast)
current_results <- topTags(current_Test,20000) 
write.table(current_results, file="Results_A_VS_F_2023-02-23.tab", sep = "\t", col.names = NA)

# B	C
current_contrast <- makeContrasts(B-C, levels=design_matrix)
current_Test <- glmQLFTest(fit, contrast=current_contrast)
current_results <- topTags(current_Test,20000) 
write.table(current_results, file="Results_B_VS_C_2023-02-23.tab", sep = "\t", col.names = NA)

# B	D
current_contrast <- makeContrasts(B-D, levels=design_matrix)
current_Test <- glmQLFTest(fit, contrast=current_contrast)
current_results <- topTags(current_Test,20000) 
write.table(current_results, file="Results_B_VS_D_2023-02-23.tab", sep = "\t", col.names = NA)

# B	E
current_contrast <- makeContrasts(B-E, levels=design_matrix)
current_Test <- glmQLFTest(fit, contrast=current_contrast)
current_results <- topTags(current_Test,20000) 
write.table(current_results, file="Results_B_VS_E_2023-02-23.tab", sep = "\t", col.names = NA)

# B	F
current_contrast <- makeContrasts(B-F, levels=design_matrix)
current_Test <- glmQLFTest(fit, contrast=current_contrast)
current_results <- topTags(current_Test,20000) 
write.table(current_results, file="Results_B_VS_F_2023-02-23.tab", sep = "\t", col.names = NA)

# C	D
current_contrast <- makeContrasts(C-D, levels=design_matrix)
current_Test <- glmQLFTest(fit, contrast=current_contrast)
current_results <- topTags(current_Test,20000) 
write.table(current_results, file="Results_C_VS_D_2023-02-23.tab", sep = "\t", col.names = NA)

# C	E
current_contrast <- makeContrasts(C-E, levels=design_matrix)
current_Test <- glmQLFTest(fit, contrast=current_contrast)
current_results <- topTags(current_Test,20000) 
write.table(current_results, file="Results_C_VS_E_2023-02-23.tab", sep = "\t", col.names = NA)

# C	F
current_contrast <- makeContrasts(C-F, levels=design_matrix)
current_Test <- glmQLFTest(fit, contrast=current_contrast)
current_results <- topTags(current_Test,20000) 
write.table(current_results, file="Results_C_VS_F_2023-02-23.tab", sep = "\t", col.names = NA)

# D	E
current_contrast <- makeContrasts(D-E, levels=design_matrix)
current_Test <- glmQLFTest(fit, contrast=current_contrast)
current_results <- topTags(current_Test,20000) 
write.table(current_results, file="Results_D_VS_E_2023-02-23.tab", sep = "\t", col.names = NA)

# D	F
current_contrast <- makeContrasts(D-F, levels=design_matrix)
current_Test <- glmQLFTest(fit, contrast=current_contrast)
current_results <- topTags(current_Test,20000) 
write.table(current_results, file="Results_D_VS_F_2023-02-23.tab", sep = "\t", col.names = NA)

# E	F
current_contrast <- makeContrasts(E-F, levels=design_matrix)
current_Test <- glmQLFTest(fit, contrast=current_contrast)
current_results <- topTags(current_Test,20000) 
write.table(current_results, file="Results_E_VS_F_2023-02-23.tab", sep = "\t", col.names = NA)


###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################

setwd("C:/Users/El3ti/Documents/Trabajo_Ciencia/Mauri_Imparalogos/Schistosoma_Protasio2012")

data <- gsub("_all_counts.tab", "" ,list.files(path=".", pattern="*_all_counts.tab"))

counts_file <- gsub(" ","", paste(data ,"_all_counts.tab"))
group_file <- gsub(" ","", paste(data ,"_grupos.tab"))

print_correlation_tab <- gsub(" ","", paste(data ,"_dif_expression.tab"))

# 1) Cargar los datos a R:
raw_counts <- read.table(counts_file, header = TRUE, sep = "\t", row.names="Term_Count")
group_list <- read.table(group_file, header = TRUE, sep = "\t")

# Groups <- c(group_list$Grupo)
# Groups

Groups <- c("A", "B", "B", "A", "A", "B", "C", "C", "C", "C")

# 2) Producción del objeto a utilizar por edgeR
edgeR_Input <- DGEList(counts=raw_counts, group=Groups)

# 3) Filtrado de la Librería:
# filterByExpr.default <- function(y, design=NULL, group=NULL, lib.size=NULL, min.count=10, min.total.count=15, ...)
keep<- filterByExpr(edgeR_Input, group=Groups)
edgeR_Input <- edgeR_Input[keep, , keep.lib.sizes=FALSE]

# 4) Normalizado
edgeR_Input<-calcNormFactors(edgeR_Input)

# 5) Diseño de matriz experimental
design_matrix <- model.matrix(~0+Groups, data=edgeR_Input$samples) # No termino de entender como funciona...
colnames(design_matrix) <- levels(edgeR_Input$samples$group)

# 6) Calculos de expresión diferencial
edgeR_Input <- estimateDisp(edgeR_Input,design_matrix)

# 7) Ajuste de la distribución binomial
fit <- glmQLFit(edgeR_Input,design_matrix)

# 8) Expresión diferencial:
# A	B
current_contrast <- makeContrasts(A-B, levels=design_matrix)
current_Test <- glmQLFTest(fit, contrast=current_contrast)
current_results <- topTags(current_Test,20000) 
write.table(current_results, file="Results_A_VS_B_2023-02-23.tab", sep = "\t", col.names = NA)

# A	C
current_contrast <- makeContrasts(A-C, levels=design_matrix)
current_Test <- glmQLFTest(fit, contrast=current_contrast)
current_results <- topTags(current_Test,20000) 
write.table(current_results, file="Results_A_VS_C_2023-02-23.tab", sep = "\t", col.names = NA)

# B	C
current_contrast <- makeContrasts(B-C, levels=design_matrix)
current_Test <- glmQLFTest(fit, contrast=current_contrast)
current_results <- topTags(current_Test,20000) 
write.table(current_results, file="Results_B_VS_C_2023-02-23.tab", sep = "\t", col.names = NA)


###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################

setwd("C:/Users/El3ti/Documents/Trabajo_Ciencia/Mauri_Imparalogos/Smansoni_SangerInstitute")

data <- gsub("_all_counts.tab", "" ,list.files(path=".", pattern="*_all_counts.tab"))

counts_file <- gsub(" ","", paste(data ,"_all_counts.tab"))
group_file <- gsub(" ","", paste(data ,"_grupos.tab"))

print_correlation_tab <- gsub(" ","", paste(data ,"_dif_expression.tab"))

# 1) Cargar los datos a R:
raw_counts <- read.table(counts_file, header = TRUE, sep = "\t", row.names="Term_Count")
group_list <- read.table(group_file, header = TRUE, sep = "\t")

# Groups <- c(group_list$Grupo)
# Groups

Groups <- c("C", "C", "C", "C", "C", "C", "A", "B", "B", "B", "A", "B", "B", "B")

# 2) Producción del objeto a utilizar por edgeR
edgeR_Input <- DGEList(counts=raw_counts, group=Groups)

# 3) Filtrado de la Librería:
# filterByExpr.default <- function(y, design=NULL, group=NULL, lib.size=NULL, min.count=10, min.total.count=15, ...)
keep<- filterByExpr(edgeR_Input, group=Groups)
edgeR_Input <- edgeR_Input[keep, , keep.lib.sizes=FALSE]

# 4) Normalizado
edgeR_Input<-calcNormFactors(edgeR_Input)

# 5) Diseño de matriz experimental
design_matrix <- model.matrix(~0+Groups, data=edgeR_Input$samples) # No termino de entender como funciona...
colnames(design_matrix) <- levels(edgeR_Input$samples$group)

# 6) Calculos de expresión diferencial
edgeR_Input <- estimateDisp(edgeR_Input,design_matrix)

# 7) Ajuste de la distribución binomial
fit <- glmQLFit(edgeR_Input,design_matrix)

# 8) Expresión diferencial:
# A	B
current_contrast <- makeContrasts(A-B, levels=design_matrix)
current_Test <- glmQLFTest(fit, contrast=current_contrast)
current_results <- topTags(current_Test,20000) 
write.table(current_results, file="Results_A_VS_B_2023-02-23.tab", sep = "\t", col.names = NA)

# A	C
current_contrast <- makeContrasts(A-C, levels=design_matrix)
current_Test <- glmQLFTest(fit, contrast=current_contrast)
current_results <- topTags(current_Test,20000) 
write.table(current_results, file="Results_A_VS_C_2023-02-23.tab", sep = "\t", col.names = NA)

# B	C
current_contrast <- makeContrasts(B-C, levels=design_matrix)
current_Test <- glmQLFTest(fit, contrast=current_contrast)
current_results <- topTags(current_Test,20000) 
write.table(current_results, file="Results_B_VS_C_2023-02-23.tab", sep = "\t", col.names = NA)
