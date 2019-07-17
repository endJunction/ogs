function printVec(msg, vec) {
    print(msg, "");
    for (i in vec)
        printf("%g ",vec[i]);
    printf("\n");
}

function sprintVec(msg, vec) {
    output = msg + " ";
    for (i in vec)
        output += sprintf("%g ",vec[i]);
    output += sprintf("\n");
    return output;
}

function printVectorIndices(name, vec, state_index) {
    for (i in vec[state_index])
    {
        printf(" s%d_%s_%d", state_index, name, i)
    }
}

function printVectorEntries(vec) {
    for (i in vec)
    {
        printf(" %g", vec[i])
    }
}

function printLine(ti) {
    printf("%g %g %g %g", timestep[ti], time[ti], dt[ti], iterations[ti]);
    printVectorEntries(g_grad[ti][0])
    printVectorEntries(g_grad[ti][1])
    printVectorEntries(g_isv[ti][0])
    printVectorEntries(g_isv[ti][1])
    printVectorEntries(g_tf[ti][0])
    printVectorEntries(g_tf[ti][1])
    printf("\n");
}
BEGIN{
    I = 0;
}
END{
    printf("timestep time dt iterations")
    printVectorIndices("grad", g_grad[1], 0)
    printVectorIndices("grad", g_grad[1], 1)
    printVectorIndices("isv", g_isv[1], 0)
    printVectorIndices("isv", g_isv[1], 1)
    printVectorIndices("tf", g_tf[1], 0)
    printVectorIndices("tf", g_tf[1], 1)
    printf("\n");
    for (i in timestep)
    {
        printLine(i);
    }
}
/Time stepping at step/ {
    #print($0);
    I++;
    split($0, line);
    timestep[I] = substr(line[7], 2);
    time[I] = line[10];
    dt[I] = line[14];
    g_grad[I][0][0] = 0;
    g_isv[I][0][0] = 0;
    g_tf[I][0][0] = 0;
    g_grad[I][1][0] = 0;
    g_isv[I][1][0] = 0;
    g_tf[I][1][0] = 0;
}
/info: \[time\] Iteration/{
    #print("Iteration", $4, "finished");
    iterations[I] = $4
}
/info: E 1997/ {
    #print("Before integration.")
}
/error: Computation of local constitutive relation failed for element 1997, integration point 1./{
    #print("Integration failed.")
}

function copyArray(input, out) {
    for (i in input)
    {
        out[i] = input[i];
    }
}

function copyState(s) {
    #print("Copy state", s, "for timestep", I);
    copyArray(grad[s], g_grad[I][s]);
    copyArray(isv[s], g_isv[I][s]);
    copyArray(tf[s], g_tf[I][s]);
    #print("Copy state finished")
}

/info: Integrate stress successful ==============/ {
    #print("Compare states.")
    #Compare s0 to s1
    #for (i in g[0])
    #{
    #    if (g[0][i] != g[1][i])
    #    {
    #        print("ERROR states not equal.")
    #        for (i in g[0])
    #            print("g0 " g[0][i] " != " g[1][i])
    #        #exit;
    #    }
    #}

    s = 2;

    #printVec("ISV", g_isv[I])
    #g_isv[I][0] = isv[0]
    #g_tf[I][0] = tf[0]

}
/info: \[time\] Assembly took/{
    #print("Integration finished");
    #print("S1 ");
    #printVec("grad", g[1])
    #printVec("isv", isv[1])
    #printVec("tf", tf[1])
    copyState(0);
    copyState(1);
}

/info: State0/{
    #print("State0");
    s = 0;}
/info: State1/{
    #print("State1");
    s = 1;}
/gradients:/{
    #print("Reading gradients");
    grad[s][0] = $3;
    grad[s][1] = $4;
    grad[s][2] = $5;
    grad[s][3] = $6;
    #printVec("grad", grad[s])
}
/internal_state_variables:/{
    #print("Reading ISV");
    isv[s][0] = $3;
    isv[s][1] = $4;
    isv[s][2] = $5;
    isv[s][3] = $6;
    isv[s][4] = $7;
    #printVec("ISV", isv[s])
}
/thermodynamic_forces:/{
    #print("Reading TF");
    tf[s][0] = $3;
    tf[s][1] = $4;
    tf[s][2] = $5;
    tf[s][3] = $6;
    #printVec("tf", tf[s])
}
