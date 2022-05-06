import os
import cobra

os.chdir("E://bioinformatics")
my_path = os.getcwd()

# setting model directory
model_dir = os.path.join(my_path, "models")
model_dict = {}
for roots, dirs, files in os.walk(model_dir):
    for each_model in files:
        model_name = os.path.splitext(each_model)[0]
        model_dict[model_name] = cobra.io.read_sbml_model(each_model)

# basic statistics
with open("basic_statistics.txt", "w+") as f1:
    for each_model in model_dict.keys():
        model = model_dict[each_model]
        f1.writelines(each_model + "\tgenes:{0}\treactions:{1}\tmetabolites:{2}".format(len(model.genes),
                                                                                        len(model.reactions),
                                                                                        len(model.metabolites)))
def extract_info(test_model):
    with open("reaction_info_{0}".format(test_model).split(".")[0] + ".txt", "w+") as f:
        gem = cobra.io.read_sbml_model(test_model)
        for each_reaction in gem.reactions:
            result_line = ''
            for a in each_reaction.reactants:
                result_line += a.id + " + "
            result_line = result_line.strip(" + ") + " <=> "
            for b in each_reaction.products:
                result_line += b.id + " + "
            result_line = result_line.strip(" + ")
            f.write(each_reaction.id + "\t" + each_reaction.gene_reaction_rule + "\t" + str(each_reaction.lower_bound)
                    + ";" + str(each_reaction.upper_bound) + "\t" + str(each_reaction.reversibility) +"\t" + str(result_line)+ '\n')
    with open("metabolite_info_{0}".format(test_model).split(".")[0] + ".txt", "w+") as outfile:
        for each_metabolite in gem.metabolites:
            outfile.write(each_metabolite.id + "\t" + each_metabolite.name + "\n")
    return 0


def generate_metabolite_in_reactions(model_file):
    reaction_metabolite_list = []
    with open(model_file, "r") as f:
        for each_line in f.readlines():
            if "speciesReference species=" in each_line:
                species = each_line.split('="M_')[1].split('"')[0]
                if species not in reaction_metabolite_list:
                    reaction_metabolite_list.append(species)
    return reaction_metabolite_list


def extract_exchange_pool(model_file):
    pool = []
    gem = cobra.io.read_sbml_model(model_file)
    for each_metabolite in gem.metabolites:
        if each_metabolite.id[-1] == 'e':
            if each_metabolite.id not in pool:
                pool.append(each_metabolite.id)
    with open("Exchange_pool.txt", "w+") as f:
        for each_exchange in pool:
            f.write(each_exchange + "\n")


def check_model(test_model):
    with open("model_check.txt", "w+") as outfile:
        gem = cobra.io.read_sbml_model(test_model)
        metabolite_list = []
        for each_metabolite in gem.metabolites:
            if each_metabolite.id not in metabolite_list:
                metabolite_list.append(each_metabolite.id)
        outfile.write("Checking reactions without gene reaction rule...\n"
                      + "------------------\n\n")  # check gene-reaction rule
        for each_reaction in gem.reactions:
            if each_reaction.gene_reaction_rule == "":
                outfile.writelines(each_reaction.id + "\n")
        outfile.write("Checking metabolites in reactions but not in metabolite list..." + "\n"
                      + "------------------\n\n")
        reaction_metabolite_list = generate_metabolite_in_reactions(test_model)
        a = [x for x in reaction_metabolite_list if x not in metabolite_list]
        b = [y for y in metabolite_list if y not in reaction_metabolite_list]
        outfile.write("Only in metabolite list\n")
        for m in b:
            outfile.write(m + "\n")
        outfile.write("Only in reaction\n")
        for n in a:
            outfile.write(n + "\n")
    return reaction_metabolite_list, metabolite_list


def get_info(test_model):
    gem = cobra.io.read_sbml_model(test_model)
    metabolite_list = []
    reaction_list = []
    for each_metabolite in gem.metabolites:
        metabolite_list.append(each_metabolite.id)
    for each_reaction in gem.reactions:
        reaction_list.append(each_reaction.id)
    return metabolite_list, reaction_list


def model_compare(model1, model2):
    metabolite1, reaction1 = get_info(model1)
    metabolite2, reaction2 = get_info(model2)

    joint_metabolite = [x for x in metabolite1 if x in metabolite2]
    joint_reaction = [y for y in reaction1 if y in reaction2]
    only1_metabolite = [z for z in metabolite1 if z not in metabolite2]
    only2_metabolite = [m for m in metabolite2 if m not in metabolite1]
    only1_reaction = [n for n in reaction1 if n not in reaction2]
    only2_reaction = [k for k in reaction2 if k not in reaction1]
    with open("Comparison_report.txt", "w+") as outfile:
        outfile.write(model1 + "only metabolites\n")
        for each_metabolite1 in only1_metabolite:
            outfile.write(each_metabolite1 + "\n")
        outfile.write(model1 + "only reactions\n")
        for each_reaction1 in only1_reaction:
            outfile.write(each_reaction1 + "\n")
        outfile.write(model2 + "only metabolites\n")
        for each_metabolite2 in only2_metabolite:
            outfile.write(each_metabolite2 + "\n")
        outfile.write(model2 + "only reactions\n")
        for each_reaction2 in only2_reaction:
            outfile.write(each_reaction2 + "\n")
        outfile.write("---------------------------------------\n")


def extract_kegg_reaction_list(model_file):
    gem = cobra.io.read_sbml_model(model_file)
    reaction_list = []
    for each_reaction in gem.reactions:
        note = each_reaction.notes
        try:
            kegg_list = note["KEGG REACTION"]
        except KeyError:
            kegg_list = []
        reaction_list += kegg_list
    with open("reaction_kegg_number_list.txt", "w+") as f:
        for each_number in list(set(reaction_list)):
            f.write(each_number + "\n")


def flux_distribution(model_file):
    gem = cobra.io.read_sbml_model(model_file)
    gem.optimize()
    reaction_number = len(gem.reactions)
    #a = flux_variability_analysis(gem, gem.reactions)
    with open("fluxes_H1S.txt", "w+") as f:
        for each_metabolite in gem.metabolites:
            f.write("--------------------------\n")
            f.write(each_metabolite.id + "\n")
            f.write("--------------------------\n")
            for each_reaction in gem.reactions:
                if each_metabolite in each_reaction.reactants:
                    f.write(each_reaction.id + ":" + str(round(-1*each_reaction.x, 4)) + "\n")
                if each_metabolite in each_reaction.products:
                    f.write(each_reaction.id + ":" + str(round(each_reaction.x, 4)) + "\n")
    return len(gem.genes), gem.optimize()


# single response plane
def response_plane(model_file):
    my_file = open("response_plane_for__ocddea.txt", "w+")
    gem = cobra.io.read_sbml_model(model_file)
    for i in [a*0.5 for a in range(0,29)]:
        for j in [a*0.5 for a in range(100)]:
            gem.reactions.get_by_id("EX_succ_e").lower_bound = -i
            gem.reactions.get_by_id("EX_succ_e").upper_bound = -i
            gem.reactions.get_by_id("EX_o2_e").upper_bound = -j
            gem.reactions.get_by_id("EX_o2_e").lower_bound = -j
            obj = gem.optimize().objective_value
            my_file.writelines(str(i) + "\t" + str(j) + "\t" + str(obj) + "\n")
    my_file.close()
    return  0
# in batch response plane
def response_plane_plus(model_file):
    exchange_reactions = []
    done_list = []
    gem = cobra.io.read_sbml_model(model_file)
    for each_reaction in gem.reactions:
        if each_reaction.id.startswith("EX_"):
            exchange_reactions.append(each_reaction.id)
    for each_metabolite in exchange_reactions:
        print(each_metabolite)
        if each_metabolite  not in done_list:
            gem = cobra.io.read_sbml_model(model_file)
            my_file = open("response_plane_{0}.txt".format(each_metabolite.split("EX_")[1]), "w+")
            gem.reactions.get_by_id(each_metabolite).lower_bound = -2
            gem.reactions.get_by_id(each_metabolite).upper_bound = -2
            for i in [a * 0.5 for a in range(0, 29)]:
                for j in [a * 0.5 for a in range(100)]:
                    gem.reactions.get_by_id("EX_succ_e").lower_bound = -i
                    gem.reactions.get_by_id("EX_succ_e").upper_bound = -i
                    gem.reactions.get_by_id("EX_o2_e").upper_bound = -j
                    gem.reactions.get_by_id("EX_o2_e").lower_bound = -j
                    my_data = gem.optimize()
                    if my_data.status == "infeasible":
                        obj = 0
                    else:
                        obj = gem.optimize().objective_value
                    my_file.writelines(str(i) + "\t" + str(j) + "\t" + str(obj) + "\n")
            my_file.close()


def extract_pathway_flux(pathway_file, model_file):
    pathway_dict = {}
    gem = cobra.io.read_sbml_model(model_file)
    gem.optimize()
    with open(pathway_file, "r") as f1:
        for each_reaction in f1.readlines():
            reaction_id = each_reaction.split("R_")[1].strip()
            flux = gem.reactions.get_by_id(reaction_id).x
            pathway_dict[reaction_id] = flux
    with open("pathway_fluxes.txt", "w+") as f2:
        for each_reaction, each_flux in pathway_dict.items():
            f2.write(each_reaction + "\t" + str(each_flux) + "\n")


def compare_flux(flux_file):
    import math
    with open("flux_compare.txt", "w+") as outfile:
        with open(flux_file, "r") as f1:
            for each_reaction in f1.readlines():
                info = each_reaction.strip().split("\t")
                flux1 = float(info[1])
                flux2 = float(info[2])
                if flux1 != 0 and flux2 != 0:
                    try:
                        logFC = math.log((flux2/flux1), 2)
                    except ValueError:
                        print(str(flux2/flux1))
                        logFC = 0
                    if logFC > 0:
                        tag = 'up'
                    elif logFC < 0:
                        tag = 'down'
                    else:
                        tag = 'constant'
                elif flux1 == 0:
                    tag = 'up'

                result_line = each_reaction.strip() + "\t" + str(logFC) + "\t" + tag
                outfile.write(result_line + "\n")


from scipy import stats
expr_dict = {}
with open("rpkm_matrix.txt", "r") as f1:
    for eachline in f1.readlines()[1:]:
        locus = eachline.strip().split("\t")[0]
        temp_list = []
        for i in eachline.strip().split("\t")[1:]:
            temp_list.append(float(i))
        expr_dict[locus] = temp_list
with open("co-expression_network.txt", "w+") as outfile:
    for j in expr_dict.keys():
        for q in expr_dict.keys():
            r, p = stats.pearsonr(expr_dict[j], expr_dict[q])
            outfile.writelines(j + "\t" + q + "\t" + str(r) + "\n")

with open("control_seq_promotors.txt", "w+") as outfile:
    with open("genome_predictions.txt") as f1:
        for eachline in f1.readlines()[1:]:
            info = eachline.strip().split("\t")
            head_line = info[0] + "_" + info[1] + "_" + info[2]
            seq = info[-1]
            outfile.writelines(">" + head_line + "\n")
            outfile.writelines(seq + "\n")

with open("reaction_bound.txt", "w+") as outfile:
    for each_reaction in gem.reactions:
        outfile.writelines(each_reaction.id + "\t" + str(each_reaction.lower_bound) + "\t" + str(each_reaction.upper_bound) + "\n")

with open("metabolite_clustering.txt_R.txt", "w+") as outfile:
    with open("metabolite_clustering.txt", "r") as f1:
        for each_line in f1.readlines():
            if each_line.startswith("Carbon sources"):
                outfile.writelines(each_line)
            else:
                info = each_line.strip().split("\t")
                cmd = ''
                for i in range(2, len(info)):
                    if float(info[i]) <= 1e-10:
                        cmd += "0\t"
                    else:
                        cmd += str(info[i])
                        cmd += "\t"
                cmd2 = cmd.strip("\t") + "\n"
                outfile.writelines(cmd2)


# statistic of flux changes for H2O2 treatments
record_dict = {}
rxn_list = []
with open("flux_alt_001mM_425.txt", "r") as f1:
        for each_line in f1.readlines():
            info = each_line.strip().split("\t")
            rxn_flux = info[1:]
            record_dict[info[0]] = rxn_flux
            if info[0].split("NF_")[1] not in rxn_list:
                rxn_list.append(info[0].split("NF_")[1])
with open("stat_flux001_reg_425.txt", "w+") as outfile:
    for each_rxn in rxn_list:
        outfile.writelines(each_rxn)
        wt_rxn = "NF_" + each_rxn
        mt_rxn = "PERTURB_NF_" + each_rxn
        for i in range(len(record_dict[wt_rxn])):
            a = float(record_dict[mt_rxn][i])- float(record_dict[wt_rxn][i])
            if a <0:
                outfile.writelines("\t" + "down")
            elif a > 0:
                outfile.writelines("\t" + "up")
            elif a==0:
                outfile.writelines("\t" +"blocked")
        outfile.writelines("\n")


# calculate connectivity of GEM metabolites
connectivity_dict = {}
gem = cobra.io.read_sbml_model("new_Azo.xml")
for each_reaction in gem.reactions:
    for each_reactant in each_reaction.reactants:
        if each_reactant not in connectivity_dict.keys():
            connectivity_dict [each_reactant] = 1
        else:
            connectivity_dict [each_reactant] += 1
    for each_product in each_reaction.products:
        if each_product not in connectivity_dict.keys():
            connectivity_dict [each_product] = 1
        else:
            connectivity_dict [each_product] += 1
with open("connectivity.txt", "w+") as outfile:
    for each_met in connectivity_dict.keys():
        outfile.writelines(each_met.id + "\t" + str(connectivity_dict[each_met]) + "\n")