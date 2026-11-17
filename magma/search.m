/* =====================================================================
   search.m

   Exhaustive parameter search for the cyclic (r,delta)-locally
   recoverable codes of the paper

       "Cyclic (r,delta) Locally Recoverable Codes from Structured
        Defining Sets".

   Construction.  N = s*n with N | q^m - 1, beta a primitive N-th root of
   unity in GF(q^m), and C the cyclic code of length N over GF(q) with
   defining set Z = Z_loc  union  Z_glob, where

       Z_loc  = union of the q-cyclotomic cosets of l*n + h  (mod N),
                0 <= l <= s-1,  0 <= h <= delta-2,
       Z_glob = union of the q-cyclotomic cosets of a + i     (mod N),
                0 <= i <= b-1.

   Search range.  Every length N <= MaxN, each at its minimal extension
   degree m = ord_N(q) <= MaxM; every factorization N = s*n with s, n >= 2;
   every 2 <= delta <= n, 0 <= a <= N-1 and b >= 1 satisfying the
   conditions of the paper, with K > r and K > r_H.  There is no
   restriction on the dimension K.

   For every code the script records K, the BCH bound d_BCH of Z, the
   locality pair (r_H, delta_H) induced by the auxiliary code B_H, the
   Singleton-type bounds B and B_H at (r, delta) and (r_H, delta_H), and
   the minimum distance d, obtained from
       - the certificate d_BCH = B_H, which gives d = B_H;
       - MinimumWeight, for N <= MaxLenD, with a limit of MaxTimeD seconds;
       - the corollary on the codes with m = 2, delta = 2, b = 1,
         n | q+1, n >= 4, gcd(a,n) = 1, for which d = 4.

   Usage.   magma -b q:=11 search.m
   Output.  <OutDir>/search_q<q>.csv  (one line per parameter set) and
            <OutDir>/search_q<q>.txt  (the same rows as LaTeX longtables,
                                       with Z_loc, Z_glob and g(x)).
   ===================================================================== */

// ---------------------------------------------------------------------
//  parameters (values given on the command line arrive as strings)
// ---------------------------------------------------------------------
if not assigned q        then q        := 3;   end if;
if not assigned MaxN     then MaxN     := 128; end if;
if not assigned MaxLenD  then MaxLenD  := 60;  end if;
if not assigned MaxTimeD then MaxTimeD := 60;  end if;
if not assigned OutDir   then OutDir   := ".";  end if;
if Type(q)        eq MonStgElt then q        := StringToInteger(q);        end if;
if Type(MaxN)     eq MonStgElt then MaxN     := StringToInteger(MaxN);     end if;
if Type(MaxLenD)  eq MonStgElt then MaxLenD  := StringToInteger(MaxLenD);  end if;
if Type(MaxTimeD) eq MonStgElt then MaxTimeD := StringToInteger(MaxTimeD); end if;

// Optional restriction of the range, used to split long runs into pieces
// that can be run in parallel and concatenated afterwards:
//   StartN[m]           : for degree m only lengths N >= StartN[m] are used;
//   SplitN, SplitFrom,  : at length N = SplitN only tuples (s,delta,b,a) with
//   SplitTo               SplitFrom < s*10^9 + delta*10^6 + b*10^3 + a <= SplitTo
//                         are written (SplitN := 0 disables this);
//   Suffix              : appended to the output file names.
if not assigned StartN    then StartN    := [0, 0, 0]; end if;
if not assigned SplitN    then SplitN    := 0;         end if;
if not assigned SplitFrom then SplitFrom := -1;        end if;
if not assigned SplitTo   then SplitTo   := 10^30;     end if;
if not assigned Suffix    then Suffix    := "";        end if;

file    := OutDir cat "/search_q" cat IntegerToString(q) cat Suffix cat ".txt";
csvfile := OutDir cat "/search_q" cat IntegerToString(q) cat Suffix cat ".csv";

MaxM         := 3;      // largest extension degree m
MinS         := 2;      // smallest number s of repair sets
MaxDimLoc    := 14;     // d_min(C|_{S_j}) is computed (if CheckRestrictions)
MaxLenLoc    := 30;     //   only for dim <= MaxDimLoc and n <= MaxLenLoc
PrintGenPoly := true;   // print g(x) in the .txt output when deg g <= MaxDegPrint
MaxDegPrint  := 18;
RequireKgtR       := true;   // keep only rows with K > r
RequireKgtRH   := true;   // and with K > r_H
CheckRestrictions := false;  // compute d_min(C|_{S_j}) for every j (slow)

// ------------------------- small helpers -----------------------------

// multiplicative order of qq modulo N  (requires gcd(qq,N) = 1)
function ModOrd(qq, N)
    e := 1;
    y := qq mod N;
    while y ne 1 do
        y := (y*qq) mod N;
        e := e + 1;
    end while;
    return e;
end function;

// q-cyclotomic coset of x modulo N
function Coset(x, qq, N)
    S := { x mod N };
    y := (x*qq) mod N;
    while y notin S do
        Include(~S, y);
        y := (y*qq) mod N;
    end while;
    return S;
end function;

// print a set of integers as a LaTeX set  \{a,b,c\}
function SetStr(S)
    L := Sort(Setseq(S));
    t := "\\{";
    for i in [1..#L] do
        t := t cat IntegerToString(L[i]);
        if i lt #L then t := t cat ","; end if;
    end for;
    return t cat "\\}";
end function;

// a set of integers as a space-separated list (CSV output)
function SpaceStr(S)
    L := Sort(Setseq(S));
    t := "";
    for i in [1..#L] do
        t cat:= IntegerToString(L[i]);
        if i lt #L then t cat:= " "; end if;
    end for;
    return t;
end function;

// LaTeX form of a polynomial: remove '*', brace the exponents
function TexPoly(f)
    s := Sprint(f);
    t := "";
    i := 1;
    digits := { "0","1","2","3","4","5","6","7","8","9" };
    while i le #s do
        c := s[i];
        if c eq "*" then
            i := i + 1;
        elif c eq "^" then
            num := "";
            j := i + 1;
            while j le #s and s[j] in digits do
                num := num cat s[j];
                j := j + 1;
            end while;
            t := t cat "^{" cat num cat "}";
            i := j;
        else
            t := t cat c;
            i := i + 1;
        end if;
    end while;
    return t;
end function;

// full BCH bound of the defining set Z:
// (longest CIRCULAR run of consecutive residues contained in Z) + 1
function BCHBound(Z, N)
    best := 0;
    cur  := 0;
    for i in [0..2*N-1] do
        if (i mod N) in Z then
            cur := cur + 1;
            if cur gt best then best := cur; end if;
        else
            cur := 0;
        end if;
    end for;
    if best gt N then best := N; end if;
    return best + 1;
end function;

// Singleton-type bound of Prakash et al. for an [N,k] code with
// (rr,dd)-locality;  returns -1 if the parameters are degenerate
function STBound(N, k, rr, dd)
    if rr lt 1 or k lt 1 then return -1; end if;
    return N - k + 1 - (Ceiling(k/rr) - 1)*(dd - 1);
end function;

// the auxiliary local code  B_H  of length n and its minimum distance.
// H is the union of the q-cyclotomic cosets of 0,...,delta-2 modulo n,
// and every restriction C|_{S_j} is contained in B_H, so
//        d_min(C|_{S_j})  >=  delta_H := d_min(B_H)  >=  delta.
// delta_H depends only on (q,n,delta), so it is cached by the caller.
// MinimumWeight(B_H) is capped by MaxT seconds.  On time-out delta_H is
// replaced by the proven lower bound max(lo, BCH bound of H) >= delta; the
// pair (n - delta_H + 1, delta_H) is then still a valid locality of C.
// Third return value: "exact", "lower", "zero" or "none".
function LocalCodeDist(q, n, delta, Fq, Pol, MaxT)
    mm := ModOrd(q, n);
    Fl := GF(q^mm);
    // no <u> here: for n | q-1 we have Fl = GF(q), and Magma returns the SAME
    // (cached) polynomial ring as P, so Pl<u> would rename x to u globally
    Pl := PolynomialRing(Fl);  u := Pl.1;
    alpha := PrimitiveElement(Fl)^((q^mm - 1) div n);
    H := {};
    for h in [0..delta-2] do
        H := H join Coset(h, q, n);
    end for;
    if #H ge n then return 0, H, "zero"; end if;      // zero code
    fl := &*[ u - alpha^h : h in H ];
    ok := true;
    for c in Eltseq(fl) do
        if c notin Fq then ok := false; break; end if;
    end for;
    if not ok then return -1, H, "none"; end if;
    f := Pol ! [ Fq ! c : c in Eltseq(fl) ];
    B := CyclicCode(n, f);
    dmin := MinimumWeight(B : MaximumTime := MaxT);
    lo, hi := MinimumWeightBounds(B);
    if lo eq hi then return lo, H, "exact"; end if;
    return Max(lo, BCHBound(H, n)), H, "lower";
end function;

// true minimum distance of the restriction C|_S  (S = 1-based positions)
// returns -1 when it is not computed
function LocalDist(C, N, S, nn, MaxDimLoc, MaxLenLoc)
    if nn gt MaxLenLoc then return -1; end if;
    T  := { 1..N } diff S;
    Cj := PunctureCode(C, T);
    if Dimension(Cj) eq 0 then return nn + 1; end if;
    if Dimension(Cj) gt MaxDimLoc then return -1; end if;
    return MinimumWeight(Cj);
end function;

// ---------------------------------------------------------------------
//  initialise the output file
// ---------------------------------------------------------------------
PrintFile(file, "% =============================================================" : Overwrite := true);
PrintFile(file, "% Cyclic (r,delta)-LRC parameter search,  q = " cat IntegerToString(q));
v1, v2, v3 := GetVersion();
PrintFile(file, "% Produced by search.m,  Magma V" cat IntegerToString(v1) cat "." cat IntegerToString(v2) cat "-" cat IntegerToString(v3));
PrintFile(file, "% Search range:  N <= " cat IntegerToString(MaxN) cat ",  m = ord_N(q) <= " cat IntegerToString(MaxM) cat
                ",  s, n >= 2,  2 <= delta <= n,  0 <= a <= N-1,  b >= 1,  K > r,  K > r_H.");
PrintFile(file, "% Columns:  t      = |Z_loc meet Z_glob|,");
PrintFile(file, "%           d_BCH  = BCH bound of the full defining set,");
PrintFile(file, "%           d      = minimum distance ('?' = not determined),");
PrintFile(file, "%           r_H,delta_H = locality pair induced by B_H,");
PrintFile(file, "%           B, G   = Singleton-type bound and gap at (r,delta),");
PrintFile(file, "%           B_H, G_H = Singleton-type bound and gap at (r_H,delta_H),");
PrintFile(file, "%           Opt    = 'A' when G_H = 0.");
if RequireKgtR then
    PrintFile(file, "% Only parameter sets with  K > r  are listed.");
end if;
PrintFile(file, "% LaTeX preamble needs:  \\usepackage{longtable}");
PrintFile(file, "% =============================================================");
PrintFile(file, "");

PrintFile(csvfile, "q,m,N,s,n,delta,r,a,b,t,K,d_BCH,d,d_lower,d_upper,d_source,r_H,delta_H,delta_H_source,B,B_H,G_H,type,Z_loc,Z_glob" : Overwrite := true);

T0     := Realtime();
Fq<w>  := GF(q);
P<x>   := PolynomialRing(Fq);

total  := 0;
totopt := 0;

// =====================================================================
//  MAIN LOOP OVER THE EXTENSION DEGREE m
// =====================================================================

for m in [1..MaxM] do

    PrintFile(file, "% -------------------- m = " cat IntegerToString(m) cat " --------------------");
    PrintFile(file, "\\begin{longtable}{ccccccccccccccc}");
    PrintFile(file, "\\caption{Construction over $\\mathbb{F}_{" cat IntegerToString(q) cat
                    "}$ with $m=" cat IntegerToString(m) cat "$: all admissible $(s,n,\\delta,a,b)$ with $k>r$.}\\\\");
    PrintFile(file, "\\hline");
    PrintFile(file, "$N$ & $s$ & $n$ & $\\delta$ & $r$ & $a$ & $b$ & $t$ & $k$ & $d_{\\mathrm{BCH}}$ & $d$ & $(r_{H},\\delta_{H})$ & $B$/$G$ & $B_{H}$/$G_{H}$ & Opt \\\\");
    PrintFile(file, "\\hline");
    PrintFile(file, "\\endhead");

    countm := 0;
    // no <W>, <y>: for m = 1, GF(q^m) and its polynomial ring are the cached
    // objects Fq and P, and naming them would rename w and x in the output
    Fbig    := GF(q^m);
    Pb      := PolynomialRing(Fbig);  y := Pb.1;
    primb   := PrimitiveElement(Fbig);

    for N in Divisors(q^m - 1) do
        if N lt 4 or N gt MaxN then continue; end if;
        if N lt StartN[m] then continue; end if;
        // treat every N at its minimal extension degree only
        if ModOrd(q, N) ne m then continue; end if;
        PrintFile(file, "% N = " cat IntegerToString(N) cat " started, elapsed "
                  cat Sprint(Round(Realtime(T0))) cat " s");

        beta := primb^((q^m - 1) div N);   // primitive N-th root of unity

        for s in Divisors(N) do
            if s lt MinS then continue; end if;
            n := N div s;
            if n lt 2 then continue; end if;

            for delta in [2..n] do

                // ---- condition (vi):  q^t*h mod n  notin {0,...,delta-2}
                okdelta := true;
                for h in [1..delta-2] do
                    for tt in [1..m-1] do
                        if ((q^tt)*h) mod n in {0..delta-2} then
                            okdelta := false;
                            break;
                        end if;
                    end for;
                    if not okdelta then break; end if;
                end for;
                if not okdelta then continue; end if;

                r := n - delta + 1;

                // ---- the auxiliary local code B_H and its distance delta_H
                deltaH, Hset, dHstatus := LocalCodeDist(q, n, delta, Fq, P, MaxTimeD);
                if deltaH le 0 then continue; end if;      // degenerate local code
                if dHstatus eq "lower" then
                    PrintFile(file, "% NOTE: d_min(B_H) not found in time at n=" cat IntegerToString(n) cat
                                    ", delta=" cat IntegerToString(delta) cat "; delta_H := lower bound "
                                    cat IntegerToString(deltaH));
                end if;
                if deltaH lt delta then
                    PrintFile(file, "% WARNING: delta_H < delta at q=" cat IntegerToString(q) cat
                                    ", n=" cat IntegerToString(n) cat ", delta=" cat IntegerToString(delta));
                end if;

                // ---- local zero set (union of q-cosets mod N)
                Zloc := {};
                for l in [0..s-1] do
                    for h in [0..delta-2] do
                        Zloc := Zloc join Coset(l*n + h, q, N);
                    end for;
                end for;

                // structure lemma:  |Z_loc| = s + m*s*(delta-2)
                if #Zloc ne s + m*s*(delta-2) then
                    PrintFile(file, "% WARNING: |Z_loc| mismatch at m=" cat IntegerToString(m) cat
                                    ", N=" cat IntegerToString(N) cat ", s=" cat IntegerToString(s) cat
                                    ", delta=" cat IntegerToString(delta) cat
                                    "  (got " cat IntegerToString(#Zloc) cat ")");
                end if;

                // ---- repair sets (1-based positions), fixed for this (s,n)
                Reps := [ { j + mu*s + 1 : mu in [0..n-1] } : j in [0..s-1] ];

                seen := {};
                for b in [1..N-1] do
                    for a in [0..N-1] do

                        I := { (a + i) mod N : i in [0..b-1] };
                        if #I ne b then continue; end if;

                        // ---- condition (iv):  (a+i)*q^t  notin  I,  1 <= t <= m-1
                        good := true;
                        for z in I do
                            for tt in [1..m-1] do
                                if ((z*q^tt) mod N) in I then
                                    good := false;
                                    break;
                                end if;
                            end for;
                            if not good then break; end if;
                        end for;
                        if not good then continue; end if;

                        Zglob := &join{ Coset(z, q, N) : z in I };
                        if #Zglob ne m*b then continue; end if;

                        Z := Zloc join Zglob;
                        k := N - #Z;
                        if k lt 1 then continue; end if;
                        if RequireKgtR and k le r then continue; end if;
                        if Z in seen then continue; end if;
                        Include(~seen, Z);
                        if N eq SplitN and
                           (s*10^9 + delta*10^6 + b*10^3 + a le SplitFrom or
                            s*10^9 + delta*10^6 + b*10^3 + a gt SplitTo) then
                            continue;
                        end if;

                        t := #(Zloc meet Zglob);

                        // ---- generator polynomial over F_{q^m}, then over F_q
                        gb := &*[ y - beta^z : z in Z ];
                        okcoef := true;
                        for c in Eltseq(gb) do
                            if c notin Fq then okcoef := false; break; end if;
                        end for;
                        if not okcoef then
                            PrintFile(file, "% WARNING: g(x) not over F_q at m=" cat IntegerToString(m) cat
                                            ", N=" cat IntegerToString(N) cat ", s=" cat IntegerToString(s) cat
                                            ", delta=" cat IntegerToString(delta) cat ", a=" cat IntegerToString(a) cat
                                            ", b=" cat IntegerToString(b));
                            continue;
                        end if;
                        g := P ! [ Fq ! c : c in Eltseq(gb) ];
                        C := CyclicCode(N, g);

                        // ---- sanity check: general dimension formula
                        if Dimension(C) ne s*n - (s + m*s*(delta-2)) - m*b + t then
                            PrintFile(file, "% WARNING: dimension mismatch at m=" cat IntegerToString(m) cat
                                            ", N=" cat IntegerToString(N) cat ", s=" cat IntegerToString(s) cat
                                            ", delta=" cat IntegerToString(delta) cat ", a=" cat IntegerToString(a) cat
                                            ", b=" cat IntegerToString(b));
                        end if;

                        dB := BCHBound(Z, N);


                        // ---- locality pair induced by B_H:
                        //      delta_H was computed once for this (n,delta)
                        dH := deltaH;
                        // optional cross-check on the actual restriction
                        if CheckRestrictions and dH ne -1 then
                            for j in [1..s] do
                                dj := LocalDist(C, N, Reps[j], n, MaxDimLoc, MaxLenLoc);
                                if dj ne -1 and dj lt dH then
                                    PrintFile(file, "% WARNING: d_min(C|_S_j) = " cat IntegerToString(dj) cat
                                                    " < delta_H = " cat IntegerToString(dH) cat
                                                    " at m=" cat IntegerToString(m) cat ", N=" cat IntegerToString(N) cat
                                                    ", s=" cat IntegerToString(s) cat ", delta=" cat IntegerToString(delta) cat
                                                    ", a=" cat IntegerToString(a) cat ", b=" cat IntegerToString(b));
                                end if;
                            end for;
                        end if;

                        // ---- the two Singleton-type bounds and gaps
                        Bnom := STBound(N, k, r, delta);
                        if dH ne -1 then
                            rH := n - dH + 1;
                            if rH lt 1 then rH := 1; end if;
                            BH := STBound(N, k, rH, dH);
                        else
                            rH := -1;
                            BH := -1;
                        end if;

                        if RequireKgtRH and rH ge 1 and k le rH then
                            continue;
                        end if;

                        // ---- minimum distance of C
                        //   dstatus = "certificate"   : d_BCH = B_H, hence d = B_H
                        //             "MinimumWeight" : computed within MaxTimeD seconds
                        //             "corollary"     : d = 4 by the corollary of the paper
                        //             "time limit"    : MinimumWeight did not finish
                        //             "not computed"  : N > MaxLenD
                        dlo := dB;  dhi := (BH ge 0) select BH else N - k + 1;
                        corollary := m eq 2 and delta eq 2 and b eq 1 and n ge 4 and
                                     (q + 1) mod n eq 0 and Gcd(a, n) eq 1;
                        if BH ge 0 and dB eq BH then
                            dtrue := BH;  dstatus := "certificate";
                        elif N le MaxLenD then
                            dtrue := MinimumWeight(C : MaximumTime := MaxTimeD);
                            lo, hi := MinimumWeightBounds(C);
                            if lo eq hi then
                                dtrue := lo;  dstatus := "MinimumWeight";
                            else
                                dtrue := -1;  dstatus := "time limit";
                                dlo := Max(dlo, lo);  dhi := Min(dhi, hi);
                            end if;
                        else
                            dtrue := -1;  dstatus := "not computed";
                        end if;
                        if dtrue lt 0 and corollary then
                            dtrue := 4;  dstatus := "corollary";
                        end if;
                        if dtrue ge 0 then dlo := dtrue; dhi := dtrue; end if;

                        if dtrue ge 0 then
                            dStr := IntegerToString(dtrue);
                            Gnom := IntegerToString(Bnom - dtrue);
                            if BH ge 0 then
                                GH   := BH - dtrue;
                                GHStr := IntegerToString(GH);
                            else
                                GH := -1; GHStr := "?";
                            end if;
                        else
                            dStr := "?"; Gnom := "?"; GH := -1; GHStr := "?";
                        end if;

                        if GH eq 0 then
                            optStr := "A";
                            totopt := totopt + 1;
                        else
                            optStr := "";
                        end if;

                        if rH ge 0 then
                            locStr := "(" cat IntegerToString(rH) cat "," cat IntegerToString(dH) cat ")";
                        else
                            locStr := "(?,?)";
                        end if;
                        if BH ge 0 then
                            BHStr := IntegerToString(BH);
                        else
                            BHStr := "?";
                        end if;

                        row := IntegerToString(N)     cat " & " cat
                               IntegerToString(s)     cat " & " cat
                               IntegerToString(n)     cat " & " cat
                               IntegerToString(delta) cat " & " cat
                               IntegerToString(r)     cat " & " cat
                               IntegerToString(a)     cat " & " cat
                               IntegerToString(b)     cat " & " cat
                               IntegerToString(t)     cat " & " cat
                               IntegerToString(k)     cat " & " cat
                               IntegerToString(dB)    cat " & $" cat
                               dStr                   cat "$ & $" cat
                               locStr                 cat "$ & $" cat
                               IntegerToString(Bnom) cat "/" cat Gnom cat "$ & $" cat
                               BHStr cat "/" cat GHStr cat "$ & " cat
                               optStr cat " \\\\";
                        PrintFile(file, row);
                        csvrow := Sprintf("%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o",
                                   q, m, N, s, n, delta, r, a, b, t, k, dB,
                                   (dtrue ge 0) select IntegerToString(dtrue) else "",
                                   dlo, dhi, dstatus, rH, dH,
                                   (dHstatus eq "lower") select "lower bound" else "exact",
                                   Bnom, BH,
                                   (dtrue ge 0 and BH ge 0) select IntegerToString(BH - dtrue) else "",
                                   (Zglob subset Zloc) select "D" else "G",
                                   SpaceStr(Zloc), SpaceStr(Zglob));
                        PrintFile(csvfile, csvrow);
                        PrintFile(file, "\\multicolumn{15}{l}{\\footnotesize $Z_{\\mathrm{loc}}=" cat SetStr(Zloc) cat
                                        "$;\\ $Z_{\\mathrm{glob}}=" cat SetStr(Zglob) cat "$} \\\\");
                        if PrintGenPoly and Degree(g) le MaxDegPrint then
                            PrintFile(file, "\\multicolumn{15}{l}{\\footnotesize $g(x)=" cat TexPoly(g) cat "$} \\\\");
                        end if;

                        countm := countm + 1;
                        total  := total + 1;
                    end for;
                end for;
            end for;
        end for;
    end for;

    PrintFile(file, "\\hline");
    PrintFile(file, "\\end{longtable}");
    if countm eq 0 then
        PrintFile(file, "% NOTE: no admissible parameters with m = " cat IntegerToString(m) cat
                        " for q = " cat IntegerToString(q) cat " in the search range.");
    end if;
    PrintFile(file, "% Rows written for m = " cat IntegerToString(m) cat " : " cat IntegerToString(countm));
    PrintFile(file, "");

end for;

PrintFile(file, "% Total rows: " cat IntegerToString(total) cat
                "   Singleton-type optimal (G_H = 0): " cat IntegerToString(totopt));
PrintFile(file, "% ----- end of search -----");

print "Search finished for q =", q;
print "Total rows:", total, "  optimal with respect to (r_H,delta_H):", totopt;
print "LaTeX tables written to:";
print file;
