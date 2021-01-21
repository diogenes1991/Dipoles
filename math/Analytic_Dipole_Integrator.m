(* Content-type: application/vnd.wolfram.mathematica *)

(*** Wolfram Notebook File ***)
(* http://www.wolfram.com/nb *)

(* CreatedBy='Mathematica 11.0' *)

(*CacheID: 234*)
(* Internal cache information:
NotebookFileLineBreakTest
NotebookFileLineBreakTest
NotebookDataPosition[       158,          7]
NotebookDataLength[      7918,        213]
NotebookOptionsPosition[      7570,        196]
NotebookOutlinePosition[      7910,        211]
CellTagsIndexPosition[      7867,        208]
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{
Cell[BoxData[
 RowBox[{
  RowBox[{"(*", 
   RowBox[{
   "This", " ", "Package", " ", "contains", " ", "Helpers", " ", "for", " ", 
    "handling", " ", "HyperGeometric", " ", "Functions"}], "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{"(*", " ", 
   RowBox[{
   "First", " ", "some", " ", "Hypergeometric", " ", "fucntion", " ", 
    "identities"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{"GaussHypergeometricExpand", "[", 
     RowBox[{"ARGS_", ",", "VAR_", ",", "ORDER_"}], "]"}], ":=", 
    RowBox[{"Module", "[", 
     RowBox[{
      RowBox[{"{", 
       RowBox[{
        RowBox[{"out", "=", "1"}], ",", 
        RowBox[{"aux", "=", "0"}]}], "}"}], ",", 
      RowBox[{
       RowBox[{"For", "[", 
        RowBox[{
         RowBox[{"III", "=", "0"}], ",", 
         RowBox[{"III", "\[LessEqual]", " ", "ORDER"}], ",", 
         RowBox[{"III", "++"}], ",", 
         RowBox[{"{", 
          RowBox[{
           RowBox[{"aux", "=", 
            RowBox[{"FullSimplify", "[", 
             RowBox[{"Sum", "[", 
              RowBox[{
               RowBox[{"FullSimplify", "[", 
                RowBox[{
                 RowBox[{"D", "[", 
                  RowBox[{
                   RowBox[{
                    FractionBox[
                    RowBox[{
                    RowBox[{"Pochhammer", "[", 
                    RowBox[{
                    RowBox[{"ARGS", "[", 
                    RowBox[{"[", "1", "]"}], "]"}], ",", "l"}], "]"}], 
                    RowBox[{"Pochhammer", "[", 
                    RowBox[{
                    RowBox[{"ARGS", "[", 
                    RowBox[{"[", "2", "]"}], "]"}], ",", "l"}], "]"}]}], 
                    RowBox[{
                    RowBox[{"Pochhammer", "[", 
                    RowBox[{
                    RowBox[{"ARGS", "[", 
                    RowBox[{"[", "3", "]"}], "]"}], ",", "l"}], "]"}], 
                    RowBox[{"l", "!"}]}]], 
                    SuperscriptBox["x", "l"]}], ",", 
                   RowBox[{"{", 
                    RowBox[{"d", ",", "III"}], "}"}]}], "]"}], "/.", 
                 RowBox[{"d", "\[Rule]", " ", "4"}]}], "]"}], ",", 
               RowBox[{"{", 
                RowBox[{"l", ",", "1", ",", "Infinity"}], "}"}]}], "]"}], 
             "]"}]}], ",", 
           RowBox[{"out", "=", 
            RowBox[{"out", "+", 
             RowBox[{"aux", 
              FractionBox[
               SuperscriptBox[
                RowBox[{"(", 
                 RowBox[{"d", "-", "4"}], ")"}], "III"], 
               RowBox[{"III", "!"}]]}]}]}]}], "}"}]}], "]"}], ";", 
       RowBox[{"Return", "[", 
        RowBox[{"Normal", "[", 
         RowBox[{"Series", "[", 
          RowBox[{
           RowBox[{
            RowBox[{
             RowBox[{
              RowBox[{
               RowBox[{"FullSimplify", "[", "out", "]"}], "/.", 
               RowBox[{
                RowBox[{"PolyLog", "[", 
                 RowBox[{"2", ",", 
                  FractionBox["1", 
                   RowBox[{"1", "-", "x_"}]]}], "]"}], "\[Rule]", " ", 
                RowBox[{
                 RowBox[{"PolyLog", "[", 
                  RowBox[{"2", ",", "x"}], "]"}], "-", 
                 FractionBox[
                  SuperscriptBox["\[Pi]", "2"], "3"], "+", 
                 RowBox[{
                  RowBox[{"Log", "[", "x", "]"}], 
                  RowBox[{"Log", "[", 
                   RowBox[{"1", "-", "x"}], "]"}]}], "-", 
                 RowBox[{
                  FractionBox["1", "2"], 
                  SuperscriptBox[
                   RowBox[{"Log", "[", 
                    RowBox[{"x", "-", "1"}], "]"}], "2"]}]}]}]}], "/.", 
              RowBox[{
               SuperscriptBox[
                RowBox[{"Log", "[", 
                 RowBox[{
                  RowBox[{"-", "1"}], "+", "x_"}], "]"}], "2"], "\[Rule]", 
               RowBox[{
                SuperscriptBox[
                 RowBox[{"Log", "[", 
                  RowBox[{"1", "-", "x"}], "]"}], "2"], "-", 
                RowBox[{"2", "I", " ", "\[Pi]", " ", 
                 RowBox[{"Log", "[", 
                  RowBox[{"1", "-", "x"}], "]"}]}], "-", 
                SuperscriptBox["\[Pi]", "2"]}]}]}], " ", "/.", 
             RowBox[{
              RowBox[{"Log", "[", 
               RowBox[{"-", "x_"}], "]"}], "\[Rule]", " ", 
              RowBox[{
               RowBox[{"Log", "[", "x", "]"}], "+", 
               RowBox[{"I", " ", "\[Pi]"}]}]}]}], "/.", 
            RowBox[{"x", "\[Rule]", " ", "VAR"}]}], " ", ",", 
           RowBox[{"{", 
            RowBox[{"d", ",", "4", ",", "ORDER"}], "}"}]}], "]"}], "]"}], 
        "]"}]}]}], "]"}]}], "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{"HyperGeometricRules", "=", 
     RowBox[{"{", 
      RowBox[{
       RowBox[{
        RowBox[{"Hypergeometric2F1", "[", 
         RowBox[{"a_", ",", "b_", ",", "c_", ",", "1"}], "]"}], "\[Rule]", 
        " ", 
        FractionBox[
         RowBox[{
          RowBox[{"Gamma", "[", "c", "]"}], 
          RowBox[{"Gamma", "[", 
           RowBox[{"c", "-", "a", "-", "b"}], "]"}]}], 
         RowBox[{
          RowBox[{"Gamma", "[", 
           RowBox[{"c", "-", "a"}], "]"}], 
          RowBox[{"Gamma", "[", 
           RowBox[{"c", "-", "b"}], "]"}]}]]}], ",", 
       RowBox[{
        RowBox[{"Hypergeometric2F1", "[", 
         RowBox[{"a_", ",", "b_", ",", "c_", ",", "z_"}], "]"}], 
        "\[RuleDelayed]", " ", 
        RowBox[{"GaussHypergeometricExpand", "[", 
         RowBox[{
          RowBox[{"{", 
           RowBox[{"a", ",", "b", ",", "c"}], "}"}], ",", "z", ",", "1"}], 
         "]"}]}]}], "}"}]}], ";"}]}]}]], "Input",
 CellChangeTimes->{{3.772365343017623*^9, 3.7723654098111486`*^9}, {
   3.772365442010744*^9, 3.772365524507777*^9}, {3.772365567928379*^9, 
   3.772365572916953*^9}, {3.772365605280472*^9, 3.772365616277701*^9}, {
   3.772365658739059*^9, 3.772365856518095*^9}, {3.772365997436513*^9, 
   3.772366050161202*^9}, {3.772366255035755*^9, 3.772366459563198*^9}, {
   3.772366552923896*^9, 3.772366557304982*^9}, {3.772366589897347*^9, 
   3.7723665910037613`*^9}, 3.7723678242286053`*^9, 3.7723678673193398`*^9, {
   3.772368014805437*^9, 3.77236802010525*^9}, {3.7723680700515203`*^9, 
   3.77236809242317*^9}, {3.772532570375445*^9, 3.7725325986653147`*^9}, {
   3.772538510047543*^9, 3.77253854643432*^9}, {3.772539581799367*^9, 
   3.7725396190649443`*^9}, 3.772886575303472*^9, {3.77288666474993*^9, 
   3.7728867027547703`*^9}, {3.774606344656671*^9, 3.774606351881509*^9}, {
   3.814794915058255*^9, 3.8147949636220703`*^9}, 3.8147950451309357`*^9}],

Cell[CellGroupData[{

Cell[BoxData[
 FractionBox[
  RowBox[{
   RowBox[{"-", "12.22990761764984"}], "-", 
   RowBox[{"(", 
    RowBox[{"-", "14.82896849316508"}], ")"}]}], 
  "7.838380179523046"]], "Input",
 CellChangeTimes->{{3.814813234692554*^9, 3.814813261549488*^9}}],

Cell[BoxData["0.33158137472140187`"], "Output",
 CellChangeTimes->{3.814813262194234*^9}]
}, Open  ]]
},
WindowSize->{1920, 1016},
WindowMargins->{{1, Automatic}, {Automatic, 11}},
FrontEndVersion->"11.0 for Linux x86 (64-bit) (September 21, 2016)",
StyleDefinitions->"Default.nb"
]
(* End of Notebook Content *)

(* Internal cache information *)
(*CellTagsOutline
CellTagsIndex->{}
*)
(*CellTagsIndex
CellTagsIndex->{}
*)
(*NotebookFileOutline
Notebook[{
Cell[558, 20, 6629, 159, 314, "Input"],
Cell[CellGroupData[{
Cell[7212, 183, 250, 7, 53, "Input"],
Cell[7465, 192, 89, 1, 66, "Output"]
}, Open  ]]
}
]
*)

