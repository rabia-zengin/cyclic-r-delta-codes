/* =====================================================================
   verify_table.m

   Independent check of the codes in results/table_best_per_length.csv
   (for every length N not dividing q-1, q+1, the Singleton-type optimal
   code of largest dimension).

   Each code is rebuilt from (q, m, N, s, n, delta, a, b):
       Z_loc  = union of the q-cyclotomic cosets of l*n + h  (mod N),
                0 <= l <= s-1, 0 <= h <= delta-2,
       Z_glob = union of the q-cyclotomic cosets of a + i     (mod N),
                0 <= i <= b-1,
       beta   = gamma^((q^m-1)/N),  gamma = PrimitiveElement(GF(q^m)),
   and the script recomputes
       K   = dim C,
       d   = d_min(C)                                (MinimumWeight),
       d_S = min_j d_min(C|_{S_j}),  S_j = { j, j+s, ..., j+(n-1)s },
       B_H = N - K + 1 - (ceil(K/r_H) - 1)(delta_H - 1),
   and checks K, d = B_H and d_S = delta_H against the table.

   Usage.  magma verify_table.m          (or, inside Magma: load "verify_table.m";)
           The result is printed and also written to verify_table_output.txt.
   ===================================================================== */

// < q, m, N, s, n, delta, a, b, K, d, r_H, delta_H >
Codes := [
    <3, 2, 8, 2, 4, 2, 1, 1, 4, 4, 3, 2>,
    <4, 2, 15, 3, 5, 2, 1, 1, 10, 4, 4, 2>,
    <5, 2, 12, 2, 6, 2, 1, 1, 8, 4, 5, 2>,
    <5, 2, 24, 4, 6, 2, 1, 1, 18, 4, 5, 2>,
    <7, 2, 12, 3, 4, 2, 1, 1, 7, 4, 3, 2>,
    <7, 2, 16, 2, 8, 2, 1, 1, 12, 4, 7, 2>,
    <7, 2, 24, 3, 8, 2, 1, 1, 19, 4, 7, 2>,
    <7, 2, 48, 6, 8, 2, 1, 1, 40, 4, 7, 2>,
    <8, 2, 63, 7, 9, 2, 1, 1, 54, 4, 8, 2>,
    <9, 3, 14, 2, 7, 3, 1, 1, 6, 5, 3, 5>,
    <9, 2, 20, 2, 10, 2, 1, 1, 16, 4, 9, 2>,
    <9, 3, 28, 4, 7, 3, 1, 1, 12, 5, 3, 5>,
    <9, 2, 40, 4, 10, 2, 1, 1, 34, 4, 9, 2>,
    <9, 3, 56, 8, 7, 3, 1, 1, 24, 5, 3, 5>,
    <9, 2, 80, 8, 10, 2, 1, 1, 70, 4, 9, 2>,
    <11, 2, 8, 2, 4, 2, 1, 1, 4, 4, 3, 2>,
    <11, 3, 14, 2, 7, 3, 1, 1, 6, 5, 3, 5>,
    <11, 2, 20, 5, 4, 2, 1, 1, 13, 4, 3, 2>,
    <11, 2, 24, 2, 12, 2, 1, 1, 20, 4, 11, 2>,
    <11, 2, 30, 5, 6, 2, 1, 1, 23, 4, 5, 2>,
    <11, 3, 35, 5, 7, 3, 1, 1, 15, 5, 3, 5>,
    <11, 2, 40, 5, 8, 2, 1, 1, 33, 4, 7, 2>,
    <11, 2, 60, 5, 12, 2, 1, 1, 53, 4, 11, 2>,
    <11, 2, 120, 10, 12, 2, 1, 1, 108, 4, 11, 2>,
    <13, 2, 21, 3, 7, 2, 1, 1, 16, 4, 6, 2>,
    <13, 2, 28, 2, 14, 2, 1, 1, 24, 4, 13, 2>,
    <13, 2, 42, 3, 14, 2, 1, 1, 37, 4, 13, 2>,
    <13, 2, 56, 4, 14, 2, 1, 1, 50, 4, 13, 2>,
    <13, 2, 84, 6, 14, 2, 1, 1, 76, 4, 13, 2>,
    <16, 2, 51, 3, 17, 2, 1, 1, 46, 4, 16, 2>,
    <16, 2, 85, 5, 17, 2, 1, 1, 78, 4, 16, 2>,
    <17, 2, 12, 2, 6, 2, 1, 1, 8, 4, 5, 2>,
    <17, 2, 24, 4, 6, 2, 1, 1, 18, 4, 5, 2>,
    <17, 2, 36, 2, 18, 2, 1, 1, 32, 4, 17, 2>,
    <17, 2, 48, 8, 6, 2, 1, 1, 38, 4, 5, 2>,
    <17, 2, 72, 4, 18, 2, 1, 1, 66, 4, 17, 2>,
    <17, 2, 96, 16, 6, 2, 1, 1, 78, 4, 5, 2>,
    <19, 2, 8, 2, 4, 2, 1, 1, 4, 4, 3, 2>,
    <19, 2, 12, 3, 4, 2, 1, 1, 7, 4, 3, 2>,
    <19, 2, 15, 3, 5, 2, 1, 1, 10, 4, 4, 2>,
    <19, 2, 24, 3, 8, 2, 1, 1, 19, 4, 7, 2>,
    <19, 2, 30, 3, 10, 2, 1, 1, 25, 4, 9, 2>,
    <19, 2, 36, 9, 4, 2, 1, 1, 25, 4, 3, 2>,
    <19, 2, 40, 2, 20, 2, 1, 1, 36, 4, 19, 2>,
    <19, 2, 45, 9, 5, 2, 1, 1, 34, 4, 4, 2>,
    <19, 2, 60, 3, 20, 2, 1, 1, 55, 4, 19, 2>,
    <19, 2, 72, 18, 4, 2, 1, 1, 52, 4, 3, 2>,
    <19, 2, 90, 9, 10, 2, 1, 1, 79, 4, 9, 2>,
    <19, 2, 120, 6, 20, 2, 1, 1, 112, 4, 19, 2>
];

function Coset(x, qq, N)
    S := { x mod N };
    y := (x*qq) mod N;
    while y notin S do Include(~S, y); y := (y*qq) mod N; end while;
    return S;
end function;

if not assigned OutFile then OutFile := "verify_table_output.txt"; end if;
nfail := 0;
head := "  q  m   N   K   d |   K   d d_S B_H | result";
print head;
PrintFile(OutFile, head : Overwrite := true);
for c in Codes do
    q, m, N, s, n, delta, a, b, K, d, rH, dH := Explode(c);
    Fq    := GF(q);
    P<x>  := PolynomialRing(Fq);
    F     := GF(q^m);
    Pb    := PolynomialRing(F);  y := Pb.1;
    beta  := PrimitiveElement(F)^((q^m - 1) div N);

    Zloc  := &join{ Coset(l*n + h, q, N) : l in [0..s-1], h in [0..delta-2] };
    Zglob := &join{ Coset(a + i, q, N) : i in [0..b-1] };
    Z     := Zloc join Zglob;

    g := P ! [ Fq ! e : e in Eltseq(&*[ y - beta^z : z in Z ]) ];
    C := CyclicCode(N, g);

    K1 := Dimension(C);
    d1 := MinimumWeight(C);

    // minimum distance of the restriction to every repair set S_j
    dl := N;
    for j in [0..s-1] do
        S  := { j + mu*s + 1 : mu in [0..n-1] };
        Cj := PunctureCode(C, { 1..N } diff S);
        dj := (Dimension(Cj) eq 0) select n + 1 else MinimumWeight(Cj);
        dl := Min(dl, dj);
    end for;

    BH := N - K1 + 1 - (Ceiling(K1/rH) - 1)*(dH - 1);
    ok := (K1 eq K) and (d1 eq d) and (dl eq dH) and (rH eq n - dH + 1) and (d1 eq BH);
    if not ok then nfail +:= 1; end if;
    line := Sprintf("%3o %2o %3o %3o %3o | %3o %3o %3o %3o | %o",
                    q, m, N, K, d, K1, d1, dl, BH, ok select "ok" else "MISMATCH");
    print line;
    PrintFile(OutFile, line);
end for;

last := Sprintf("%o codes checked, %o mismatches.", #Codes, nfail);
print "";
print last;
PrintFile(OutFile, "");
PrintFile(OutFile, last);
