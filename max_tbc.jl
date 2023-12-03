include("src/APXTBC.jl")


datasets = [
"16_brain_100206_90.txt",
"17_brain_100206_70.txt",
"01_hypertext.txt",
"02_highschool.txt",
"03_hospital_ward.txt",
"04_college_msg.txt",
"05_wiki_elections.txt",
"06_highschool.txt",
"07_digg_reply.txt",
"08_infectious.txt",
"09_primary_school.txt",
"10_facebook_wall.txt",
"11_slashdot_reply.txt",
"12_highschool.txt",
"13_topology.txt",
"14_SMS.txt",
"18_venice.txt",
"19_bordeaux.txt",
"21_mathoverflow.txt",
"20_askubuntu.txt"

]

path_opt = ["sfm"]
for p in path_opt
    get_max_temporal_bc(p,datasets)
end
