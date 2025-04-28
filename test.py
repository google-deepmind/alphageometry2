# Copyright 2025 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ==============================================================================

"""Test AlphaGeometry's logic core on IMO problems."""

from ddar import DDAR
from parse import AGProblem


problems_without_aux = {
    # TODO(miroslavo): Add here.
}

problems_with_aux = {
    "2001_p5a": (
        "a@0.3298517509465857_1.4439692260026906 = ;"
        " b@0.7379367796811421_0.6243724929135032 = ;"
        " c@-0.44500078968137113_0.27466310281841355 = ;"
        " p@0.2707582474456216_0.4862614684562726 = ;"
        " q@-0.00025251512413585993_0.9458190369448746 = ;"
        " r@-0.1758973264001069_0.6807588576698682 = ;"
        " t@0.9550742881108957_0.18827418003981203 = ;"
        " u@-0.4450007896813813_0.2746631028183983 = ;"
        " x@0.4792133626001627_0.6403365205958365 = ;"
        " o@0.5866446683502685_0.276936133597599 = aconst a b a c 60, eqangle a"
        " b a p a p a c, coll p b c, eqangle b a b q b q b c, coll q a c,"
        " distseq a b b p a q q b 1 1 -1 -1, coll r a c, perp r b a p, coll t a"
        " b, cong b t b p, coll u a c, perp u t a p, eqangle t b t x t x t p,"
        " coll x b r, cong o b o t, cong o t o p ? acompute b a b c"
    ),
    "2003_p4b": (
        "a@-2.0_0.0 = ; b@2.0_0.0 = ; c@1.9109057743109794_3.051295590654324 ="
        " ; o@0.0_1.4685508816796977 = ;"
        " d@-0.9732806727898079_3.7509545322323428 = ;"
        " p@1.8880370834465432_3.834501629643334 = ;"
        " q@0.4573782053283676_1.917250814821666 = ; r@-0.9732806727898079_0.0"
        " = ; x@0.21813684153819998_1.7305942803610102 = ;"
        " x'@0.21813684153820143_1.7305942803610088 = ;"
        " xb@0.6401728961839753_1.3207013102933183 = ;"
        " xd@-0.06405157724976585_2.209118607896344 = cong o a o b, cong o b o"
        " c, cong o a o d, coll p b c, perp p d b c, coll q c a, perp q d c a,"
        " coll r a b, perp r d a b, cong p q q r, eqangle b a b x b x b c,"
        " eqangle d a d x d x d c, coll x' a c, eqangle b a b x' b x' b c, coll"
        " xb b x', cong c xb c x', coll xd d x', cong c xd c x' ? coll x a c"
    ),
    "2005_p1": (
        "a@0.0_3.4641016151377544 = ; b@-2.0_0.0 = ; c@2.0_0.0 = ;"
        " a1@-0.4350517329096339_0.0 = ; a2@0.9469741497719774_0.0 = ;"
        " b1@1.217525866454817_1.3552849549086916 = ;"
        " b2@0.5265129251140116_2.552154477998579 = ;"
        " c1@-0.7824741335451829_2.108816660229063 = ;"
        " c2@-1.473487074885988_0.9119471371391762 = ;"
        " x@-0.0_1.1547005383792515 = ;"
        " p@-0.16450001622679433_1.3552849549086916 = cong a b b c, cong b c c"
        " a, coll a1 b c, coll a2 b c, coll b1 c a, cong a1 a2 a2 b1, coll b2 c"
        " a, cong a2 b1 b1 b2, coll c1 a b, cong b1 b2 b2 c1, coll c2 a b, cong"
        " b2 c1 c1 c2, cong c1 c2 c2 a1, coll x a1 b2, coll x b1 c2, para p a1"
        " a2 b1, para p b1 b c ? coll x c1 a2"
    ),
    "2008_p1b": (
        "a@0.0_0.0 = ; b@1.0_0.0 = ; c@0.290854532548588_0.7063937637795393 = ;"
        " h@0.29085453254858795_0.2919875344608819 = ;"
        " d@0.645427266274294_0.35319688188976966 = ;"
        " e@0.145427266274294_0.35319688188976966 = ; f@0.5_0.0 = ;"
        " a1@0.3905039838452928_0.6071309831641247 = ;"
        " a2@0.9003505487032953_0.0992627806154146 = ;"
        " b1@0.2055009708273334_0.49909692990195786 = ;"
        " b2@0.08535356172125463_0.20729683387758152 = ;"
        " c1@0.14083632305588817_-0.0 = ; c2@0.8591636769441118_-0.0 = ;"
        " x@0.20720154868679602_0.20800868670467504 = ;"
        " y@0.854104386031649_0.060071878854946054 = ;"
        " z@0.29085453254858795_0.4144062293186577 = ;"
        " t@0.09067630642791989_0.14097780649808872 = perp h a b c, perp h b c"
        " a, perp h c a b, coll d b c, cong d b d c, coll e a c, cong e a e c,"
        " coll f a b, cong f a f b, cong d a1 d h, coll a1 b c, cong d a2 d h,"
        " coll a2 b c, cong e b1 e h, coll b1 c a, cong e b2 e h, coll b2 c a,"
        " cong f c1 f h, coll c1 a b, cong f c2 f h, coll c2 a b, cong e x e"
        " b1, cong f x f c1, cong f y f c1, cong d y d a1, cong d z d a1, cong"
        " e z e b1, coll t a a1, cyclic a1 a2 b1 t ? cyclic c1 c2 b1 a1"
    ),
    "2008_p6": (
        "x@4.96_-0.13 = ; y@-1.006896832888816_-1.253488108068277 = ;"
        " z@-2.840284723857512_-4.911776273400683 = ;"
        " o@2.8799999999999986_-5.489999999999999 = ;"
        " w@6.909004923003877_-1.3884003936987552 = ;"
        " a@-2.080746196980964_2.602229867485152 = ;"
        " b@-1.6356942531718726_7.005065168739374 = ;"
        " c@3.015094469495374_2.4365912907211023 = ;"
        " d@1.627271672827975_1.1632975597981 = ;"
        " i1@-0.5336859138585403_3.955774130812619 = ;"
        " i2@1.4937014337487793_1.8726952212475092 = ;"
        " f1@-0.5792688026773871_2.5534248527852306 = ;"
        " f2@1.5136170751917972_2.4853963054210246 = ;"
        " q@0.6999287191540016_4.624247192763765 = ;"
        " t@2.032680584622158_2.16475810188334 = ;"
        " p@-1.2353352672275948_2.7407229046275665 = ;"
        " s@1.1871434839892263_1.3418266172024895 = ;"
        " k@3.0667847835076225_0.25639986814788507 = ;"
        " p'@-0.4881030250396932_5.358123408840009 = ;"
        " ci@-1.929662061825108_4.096883633134149 = ;"
        " ai@0.44955094305798904_4.956726967133048 = ;"
        " e2@1.388524878059561_-1.363042371478294 = ;"
        " ce@-2.442424147293411_-0.9757984266772713 = ;"
        " ae@4.086801839440497_1.3838521928104057 = ;"
        " c2@1.2719246196873397_1.3011934311661078 = ;"
        " a2@1.9081358549908851_1.42098362984077 = ;"
        " q'@1.4737857923057616_1.2599941370739938 = ;"
        " e1@-0.2929579785797208_11.361725245075798 = ;"
        " ce2@-3.4812611312229262_3.145713274802929 = ;"
        " ae2@5.665025399448304_4.867838600241024 = ;"
        " k'@3.0667847835076225_0.25639986814788507 = cong o x o y, cong o y o"
        " z, cong o w o x, perp a z o z, perp a x o x, perp b z o z, perp b w o"
        " w, perp c y o y, perp c w o w, perp d x o x, perp d y o y, eqangle a"
        " b a i1 a i1 a c, eqangle c a c i1 c i1 c b, eqangle b c b i1 b i1 b"
        " a, eqangle a c a i2 a i2 a d, eqangle d a d i2 d i2 d c, eqangle c d"
        " c i2 c i2 c a, perp f1 i1 a c, coll f1 a c, perp f2 i2 a c, coll f2 a"
        " c, cong i1 q i1 f1, cong i2 t i2 f2, perp q i1 q t, perp t i2 t q,"
        " cong i1 p i1 f1, cong i2 s i2 f2, perp p i1 p s, perp s i2 s p, coll"
        " k q t, coll k p s, coll p' f1 i1, cong i1 p' i1 f1, perp ci i1 a b,"
        " coll ci a b, perp ai i1 b c, coll ai b c, eqangle b c b e2 b e2 b a,"
        " eqangle a b a e2 a e2 a c, eqangle c a c e2 c e2 c b, perp ce e2 a b,"
        " coll ce a b, perp ae e2 b c, coll ae b c, perp c2 i2 a d, coll c2 a"
        " d, perp a2 i2 c d, coll a2 c d, coll q' i2 f2, cong i2 q' i2 f2,"
        " eqangle d a d e1 d e1 d c, eqangle c d c e1 c e1 c a, eqangle a c a"
        " e1 a e1 a d, perp ce2 e1 a d, coll ce2 a d, perp ae2 e1 c d, coll ae2"
        " c d, cong o k' o x, para k' o i1 f1 ? cong o k o x"
    ),
    "2009_p4a": (
        "a@-0.13106393182449125_0.991373918244123 = ;"
        " b@0.9235495324287674_-0.383479153475916 = ;"
        " c@-0.7944148862784848_-0.6093746688323302 = ;"
        " d@0.06456732307514113_-0.4964269111541231 = ;"
        " e@-0.46273940905148797_0.1909996247058966 = ;"
        " k@-0.29118373555788163_-0.22335934974252558 = ;"
        " i@-0.0006430952247363558_-0.0004933013547075463 = ;"
        " p@0.023225574464428453_-0.1820176011318129 = ;"
        " r@-0.584141059528896_-0.1019576992651177 = ;"
        " f@-0.46273940905148797_0.1909996247058966 = ;"
        " s@-0.16978208508047374_0.06959797422848871 = cong a b a c, coll d b"
        " c, eqangle a b a d a d a c, coll e a c, eqangle b a b e b e b c,"
        " eqangle d a d k d k d c, eqangle c d c k c k c a, eqangle e b e k d k"
        " d a, coll i a d, coll i b e, perp p k a d, coll p a d, perp r k a c,"
        " coll r a c, perp f i a c, coll f a c, perp s k i f, coll s i f ?"
        " acompute a b a c"
    ),
    "2009_p4b": (
        "a@-0.13106393182449125_0.991373918244123 = ;"
        " b@0.9235495324287674_-0.383479153475916 = ;"
        " c@-1.5059170035445306_-0.06323954600913563 = ;"
        " d@-0.2911837355578814_-0.22335934974252583 = ;"
        " e@-0.7005467204012409_0.5545387182891501 = ;"
        " k@-0.6000728677070244_0.1793257918291189 = ;"
        " i@-0.22485994124699302_0.2797996445233356 = ;"
        " p@-0.24428573084663047_0.13242778711786796 = ;"
        " r@-0.8184904676845108_0.4640671861174937 = ;"
        " f@-0.5337490733961359_0.6824847860949803 = ;"
        " s@-0.31533147341864937_0.39774339180660556 = cong a b a c, coll d b"
        " c, eqangle a b a d a d a c, coll e a c, eqangle b a b e b e b c,"
        " eqangle d a d k d k d c, eqangle c d c k c k c a, eqangle e b e k d k"
        " d a, coll i a d, coll i b e, perp p k a d, coll p a d, perp r k a c,"
        " coll r a c, perp f i a c, coll f a c, perp s k i f, coll s i f ?"
        " acompute a b a c"
    ),
    "2011_p6": (
        "a@0.0_0.0 = ; b@1.0283476899183277_0.034528533896763136 = ;"
        " c@0.40058612850850084_0.78168033091214 = ;"
        " o@0.5069062850742152_0.23371071212254715 = ;"
        " p@0.9656829052769494_0.5516709164939837 = ;"
        " q@0.8162546035654066_0.767277178785135 = ;"
        " pa@0.529747243327681_0.18539507681495893 = ;"
        " pb@-0.11614279998371435_1.1060719342164795 = ;"
        " pc@1.000512815552948_-0.4856525296585569 = ;"
        " qa@0.3431250264084389_0.36975078148314655 = ;"
        " qb@0.14614247983733242_1.110688175192835 = ;"
        " qc@0.8658834578006742_-0.7107965140835621 = ;"
        " a1@1.9743170244664063_1.1428641869408247 = ;"
        " b1@1.0779509857836376_-0.35615077399839645 = ;"
        " c1@-0.39724046699762716_1.1011245937873375 = ;"
        " o1@0.7936534995773209_0.8313581049087352 = ;"
        " x@0.2654452205404853_-0.2695498894774664 = ;"
        " m@0.2654452205404853_-0.2695498894774664 = cong o a o b, cong o b o"
        " c, cong o p o a, perp q p o p, cong b p b pa, cong c p c pa, perp b c"
        " p pa, cong c p c pb, cong a p a pb, perp c a p pb, cong a p a pc,"
        " cong b p b pc, perp a b p pc, cong b q b qa, cong c q c qa, perp b c"
        " q qa, cong c q c qb, cong a q a qb, perp c a q qb, cong a q a qc,"
        " cong b q b qc, perp a b q qc, coll a1 pb qb, coll a1 pc qc, coll b1"
        " pa qa, coll b1 pc qc, coll c1 pa qa, coll c1 pb qb, cong o1 a1 o1 b1,"
        " cong o1 b1 o1 c1, cong o x o a, coll x o o1, cyclic b1 b pa m, cyclic"
        " b1 a1 c1 m ? cong o1 a1 o1 x"
    ),
    "2013_p3": (
        "a@0.07377274994178976_-0.997275078083287 = ;"
        " b@-0.9291456762431524_0.36971382489522236 = ;"
        " c@0.9187738953521658_-0.37732328693639766 = ;"
        " ia@0.5823388125357152_2.313222520403771 = ;"
        " a1@-0.3052939208820899_0.11751647783553433 = ;"
        " ib@1.0143630192653978_-1.1417706638836451 = ;"
        " b1@0.616315677503773_-0.5992277363216619 = ;"
        " ic@-1.7717719119880044_-0.7137583697528488 = ;"
        " c1@-0.7072412268578879_0.06725560704682916 = ;"
        " o@-0.37870444636130307_-0.9277645168182471 = ;"
        " b0@-0.5947165497261434_0.7997320753254628 = ;"
        " c0@0.7983509159005565_0.5857259282600635 = eqangle b c b ia b ia b a,"
        " eqangle c b c ia c ia c a, coll a1 b c, perp ia a1 b c, eqangle a c a"
        " ib a ib a b, eqangle c a c ib c ib c b, coll b1 a c, perp ib b1 a c,"
        " eqangle b a b ic b ic b c, eqangle a b a ic a ic a c, coll c1 a b,"
        " perp ic c1 a b, cong o a1 o b1, cong o a1 o c1, cyclic o a b c, coll"
        " b0 ia ic, cyclic a b c b0, coll c0 ia ib, cyclic a b c c0 ? perp a"
        " b a c"
    ),
    "2019_p2": (
        "a@0.0_0.0 = ; b@1.0_0.0 = ; c@-0.7858234607409786_1.453721281271722 ="
        " ; a1@-0.294791634188577_1.0540046061728465 = ;"
        " b1@-0.25967235866767224_0.4803766403172407 = ;"
        " p@-0.22206498335262229_0.7939761111865026 = ;"
        " q@-1.082011648284123_0.7939761111865026 = ;"
        " p1@-0.1109779712303118_1.7203057986524242 = ;"
        " q1@0.003897850865149388_1.152665432722081 = ;"
        " o@0.5_1.2095329349302175 = ; a2@-0.554946252758179_1.9841672515418574"
        " = ; b2@-0.6784220374653078_0.6400670252420094 = coll b c a1, coll a c"
        " b1, coll a a1 p, coll b b1 q, para a b p q, coll b1 p p1, eqangle a b"
        " a c p1 p p1 c, coll a1 q q1, eqangle b a b c q1 q q1 c, cong o a o b,"
        " cong o b o c, cong o a2 o a, coll a2 a a1, cong o b2 o b, coll b2 b"
        " b1 ? cyclic p q p1 q1"
    ),
    "2021_p3": (
        "a@0.35158228874560216_0.6253011491766167 = ;"
        " b@-0.32787335987583544_-0.11434443729989269 = ;"
        " c@0.593611241171965_-0.13630297414077203 = ;"
        " d@0.25541529417090475_0.19047220899818684 = ;"
        " e@0.4472680874521529_0.32420205294349663 = ;"
        " f@0.18152091649019214_0.4401748003500804 = ;"
        " x@0.1647645897963786_1.2131693707344662 = ;"
        " o1@0.7063662227386552_0.3187883084445305 = ;"
        " o2@0.0513262541690787_0.687748158526416 = ;"
        " y@1.5549937258179136_-0.15921225844489367 = ;"
        " da@0.24023420886499053_-0.44659830430440767 = ;"
        " db@0.6811243375628278_0.3257576028873138 = ;"
        " dc@-0.07354195379055742_0.49265999240284974 = ;"
        " q@0.22915675122419496_0.07174151089699222 = ;"
        " m@0.49881279704621323_0.5293204195487172 = ;"
        " k@0.5823757264495235_0.7709477031088625 = ;"
        " y'@1.5549937258179136_-0.15921225844489367 = ;"
        " e'@0.4472680874521529_0.32420205294349663 = ;"
        " p@0.4052374543908715_0.4564618636983748 = ;"
        " m'@0.49881279704621323_0.5293204195487172 = eqangle a b a d a d a c,"
        " eqangle d e d a c d c b, coll e a c, eqangle d f d a b d b c, coll f"
        " a b, cong x b x c, eqangle b x b c c b c x, coll x a c, cong o1 a o1"
        " d, cong o1 d o1 c, cong o2 e o2 x, cong o2 x o2 d, coll y e f, coll y"
        " b c, cong b d b da, cong c d c da, perp b c d da, cong c d c db, cong"
        " a d a db, perp c a d db, cong a d a dc, cong b d b dc, perp a b d dc,"
        " cong q da q db, cong q db q dc, coll m a y, cyclic a b c m, cong y k"
        " y d, cyclic a d c k, coll y' b c, eqangle d y' d c b d b c, coll e' f"
        " y', cyclic f b c e', coll p a c, coll p d k, coll m' p b, cyclic b d"
        " k m' ? coll o1 o2 y"
    ),
}

explanation_without_aux = """
We run the logical core DDAR on some easier IMO problems that can be solved by DDAR alone.
"""

explanation_with_aux = """
We run the logical core DDAR on some challenging IMO problems, with manually provided
auxiliary points. This is not a full AlphaGeometry system, only a test of the logical core.
"""


def print_problem_and_solve(problems_dict: dict[str, str]) -> None:
  """Prints problem ID and its proving status."""
  for name, pstring in problems_dict.items():
    print("Problem:", name)
    problem = AGProblem.parse(pstring)
    ddar = DDAR(problem.points)
    for pred in problem.preds:
      ddar.force_pred(pred)
    ddar.deduction_closure()
    if ddar.check_pred(problem.goal):
      print(" Proven :-)")
    else:
      print()
      print()
      print("!!! Problem not solved, missing an auxiliary point?")
    print()


if __name__ == "__main__":
  # Easier problems that can be solved without auxiliary points.
  print(explanation_without_aux)
  print_problem_and_solve(problems_without_aux)

  # Harder problems supplemented by auxiliary points found by a language model.
  print(explanation_with_aux)
  print_problem_and_solve(problems_with_aux)
