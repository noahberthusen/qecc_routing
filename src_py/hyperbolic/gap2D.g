LoadPackage("LINS");

r := 5;;
s := 4;;

F := FreeGroup("a","b");;
G := F / [ F.1^r, F.2^s, (F.1*F.2)^2 ];;

L := LowIndexNormalSubgroupsSearchForAll(G, 1000);;

for H in List(L) do
    Gh := G / Grp(H);
    GhFp := Image(IsomorphismFpGroupByGenerators(Gh, GeneratorsOfGroup(Gh)));;

    rhos := Subgroup(GhFp, [GhFp.1]);;
    rhos_left := LeftCosets(GhFp, rhos);;
    sigmas := Subgroup(GhFp, [GhFp.2]);;
    sigmas_left := LeftCosets(GhFp, sigmas);;
    taus := Subgroup(GhFp, [GhFp.2*GhFp.1]);;
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