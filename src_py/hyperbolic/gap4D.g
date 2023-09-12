LoadPackage("LINS");

p := 5;;
q := 3;;
r := 3;;
s := 5;;

F := FreeGroup("a","b","c","d","e");;
G := F / [ F.1^2, (F.1*F.2)^p, (F.1*F.3)^2, (F.1*F.4)^2, (F.1*F.5)^2, F.2^2, (F.2*F.3)^q, (F.2*F.4)^2, (F.2*F.5)^2, F.3^2, (F.3*F.4)^r, (F.3*F.5)^2, F.4^2, (F.4*F.5)^s, F.5^2];;
# G := F / [ F.1^2, (F.1*F.2)^p, F.2^2, (F.2*F.3)^q, F.3^2, (F.3*F.4)^r, F.4^2, (F.4*F.5)^s, F.5^2];;

L := LowIndexNormalSubgroupsSearchForAll(G, 100000);;

for H in List(L) do
    Gh := G / Grp(H);
    GhFp := Image(IsomorphismFpGroupByGenerators(Gh, GeneratorsOfGroup(Gh)));;

    rhos := Subgroup(GhFp, [GhFp.1, GhFp.4, GhFp.3, GhFp.5]);;
    rhos_left := LeftCosets(GhFp, rhos);;
    sigmas := Subgroup(GhFp, [GhFp.1, GhFp.2, GhFp.3, GhFp.5]);;
    sigmas_left := LeftCosets(GhFp, sigmas);;
    taus := Subgroup(GhFp, [GhFp.1, GhFp.2, GhFp.4, GhFp.5]);;
    taus_left := LeftCosets(GhFp, taus);;

    str_name := Concatenation("../", String(r), "_", String(s), "_", String(Size(taus_left)), ".txt");;
    output := OutputTextFile(str_name, false);;
    SetPrintFormattingStatus(output, false);


    PrintTo(output, "r,");
    PrintTo(output, r);
    PrintTo(output, "\n");

    PrintTo(output, "s,");
    PrintTo(output, s);
    PrintTo(output, "\n");

    PrintTo(output, "n,");
    PrintTo(output, Size(taus_left));
    PrintTo(output, "\n");

    PrintTo(output, "m1,");
    PrintTo(output, Size(rhos_left));
    PrintTo(output, "\n");

    PrintTo(output, "m2,");
    PrintTo(output, Size(sigmas_left));
    PrintTo(output, "\n");


    for rlc in rhos_left do
        i := 0;
        for tlc in taus_left do
            if Size(Intersection(tlc, rlc)) = 1 then
                PrintTo(output, i);
                PrintTo(output, ",");
            fi;
            i := i + 1;
        od;
        PrintTo(output, "\n");
    od;

    for slc in sigmas_left do
        i := 0;
        for tlc in taus_left do
            if Size(Intersection(tlc, slc)) = 1 then
                PrintTo(output, i);
                PrintTo(output, ",");
            fi;
            i := i + 1;
        od;
        PrintTo(output, "\n");
    od;
od;