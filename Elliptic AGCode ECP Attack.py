# https://www.desmos.com/calculator/rz1syoxfie
def DecodeECP(y, A, B, C):
    code = codes.LinearCode(C)
    l, m = A.nrows(), B.nrows()
    k, n = C.dimensions()
    temp = []
    for i in range(m):
        s = []
        for j in range(l): s.append(y.dot_product(A[j].pairwise_product(B[i])))
        temp.append(s)
    T = matrix(temp)
    M1 = T.right_kernel().basis()
    nonzeroA = vector(F, M1[0]) * A
    J = [el for el in range(n) if nonzeroA[el] == 0]
    H = code.parity_check_matrix()
    Hj = H.delete_columns([i for i in range(H.ncols()) if i not in J])
    eq = H * matrix(y).transpose()
    res = vector(Hj.solve_right(eq))
    e = [0] * n
    for i in range(len(J)):
        e[J[i]] = res[i]
    c = vector([y[j] - e[j] for j in range(n)])
    return c, e if c in code else "Error", _

def RecoverCode(subcode):
    s_squared = SquareCode(subcode)
    sDual_squared = codes.LinearCode(s_squared).parity_check_matrix()
    T = MultiplyCode(subcode, sDual_squared)
    return codes.LinearCode(T).parity_check_matrix()

def CheckECP(A, B, C, t, j):
    Acode = codes.LinearCode(A)
    Bdcode = codes.LinearCode(B).dual_code()
    C_code = codes.LinearCode(C)
    d1 = MultiplyCode(A, B) == C_code.parity_check_matrix().delete_columns([j]).echelon_form()
    d2 = Acode.dimension() > t
    #d3 = Bdcode.minimum_distance(algorithm="guava") > t
    #d4 = Acode.minimum_distance(algorithm="guava") + C_code.minimum_distance(algorithm="guava") > n
    return d1 and d2 #and d3 and d4
    
def SquareCode(code):
    squared = []
    k, n = code.dimensions()
    for i in range(k):
        for j in range(i, k):
                squared.append(code[i].pairwise_product(code[j]))
    M = matrix(F, squared)
    return codes.LinearCode(M).generator_matrix()

def MultiplyCode(A, B):
    multiplied = []
    for i in range(A.nrows()):
        for j in range(B.nrows()):
            multiplied.append(A[i].pairwise_product(B[j]))
    return codes.LinearCode(matrix(F, multiplied)).generator_matrix().echelon_form()
            
def AttackECP(C):
    C = C.echelon_form()
    k, n = C.dimensions()
    H = codes.LinearCode(C).parity_check_matrix()
    H_s = SquareCode(H)
    k1, k2 = H.nrows(), H_s.nrows()
    deg, g = k2 - k1, k2 - 2 * k1 + 1
    if g < 0: return "Error", "g", _, _
    print("deg(G) =", deg, "\ng =", g)
    t = get_t(deg, g)
    if t < 0: return "Error", "t", _, _
    print("t = {0}\n".format(t))
    for i in range(n):
        if len(H.nonzero_positions_in_column(i)) == 1:
            j = i
            break
    FILTRATION = []
    FILTRATION.append(H.delete_columns([j]).echelon_form())
    FILTRATION.append(H.delete_columns([j]).delete_rows([j]).echelon_form())
    for i in range(2, t + g + 1):
        Hi = codes.LinearCode(FILTRATION[i - 1]).parity_check_matrix()
        Bi_s = SquareCode(FILTRATION[i - 1])
        Hi_s = codes.LinearCode(Bi_s).parity_check_matrix()
        if Hi_s.nrows() == 0: 
            return "Error", "Hi squared", _, _
        temp = MultiplyCode(FILTRATION[i - 2], Hi_s)
        Bi = matrix(F, temp.stack(Hi).echelon_form().right_kernel().basis())
        FILTRATION.append(Bi)
    B0_dual = codes.LinearCode(FILTRATION[0]).parity_check_matrix()
    A = codes.LinearCode(MultiplyCode(FILTRATION[-1], B0_dual)).parity_check_matrix()
    return (A, FILTRATION[-1], FILTRATION[0], j) if CheckECP(A, FILTRATION[-1], C, t, j) else ("Error", "ECP", _, _)

def TraceCodeDelsarte(code):
    C = codes.LinearCode(code).dual_code()
    tr = codes.SubfieldSubcode(C, GF(p)).dual_code()
    return tr.generator_matrix()

def TraceCodeMatrix(code):
    k, n = code.dimensions()
    tr = []
    for i in range(k):
        tr.append([code[i][e].trace() for e in range(n)])
        for j in range(i, k):
            tr.append([(code[i][e] + a * code[j][e]).trace() for e in range(n)])
    TrM = matrix(tr)
    return codes.LinearCode(TrM).generator_matrix()

def dim_Tr(m, p, k):
    #if k < (p ** (m / 2) - 3):
    #    return m * (k - floor(k / p)) + 1
    return m * (k - floor(k / p)) + 1

def encode_Tr(u, code_matrix):
    return vector([e.trace() for e in u * code_matrix])

def get_t(deg, g):
    return floor((get_d(deg, g) - g - 1) / 2)

def get_d(deg, g):
    return deg - 2 * g + 2

###################### Задаем проективное пространство с координатами над расширением поля
p, m = 5, 2
q = p ** m
F.<a> = GF(q)
P.<x,y,z> = ProjectiveSpace(F, 2)

###################### Задаем кривую
f = EllipticCurve(F, [2, 1])
EC = Curve(f)
g = EC.genus()

###################### Строим дивизоры D и G
FF = EC.function_field()
D = FF.places()
###################### Задаем размерность кода
Q, deg = EC([0, 1, 0]).place(), 6
######################
D.remove(Q)
G = deg * Q

###################### Строим АГ-код, выводим базис пространства Римана-Роха
C = codes.EvaluationAGCode(D, G)
RRBasis = C.basis_functions()

n, k, d_, t = len(D), deg, get_d(deg, g), get_t(deg, g)
print("E:", f)
print("#E =", n, "\ng =", g, "\n")
print("Basis of Riemann-Roch space:", RRBasis)
print("G = {0} * Q = {0} * {1}".format(deg, Q), "\n")

print("C(D, G):")
print(C)
CMatrix = C.generator_matrix()
print(CMatrix)
print()

###################### Строим дуальный код
C_Dual, k_dual = C.dual_code(), n - k
print("C(D, G) Dual:")
print(C_Dual)
CDualMatrix = C_Dual.generator_matrix()
print(CDualMatrix)
print("\nd* =", d_)
print("t =", t)
print()
print("--------------------------------------------------------------------------------------------------------\n")

###################### Строим трейс-код
print("dim(Tr(C)) =", dim_Tr(m, p, k))
TrCode = TraceCodeDelsarte(CMatrix)#codes.LinearCode(TraceCodeDelsarte(CMatrix)).parity_check_matrix()
print("Tr(C):")
print(TrCode)
k_tr = TrCode.nrows()
print("dim(Tr(C)) =", k_tr, "\n")

###################### Атака
A, B, PuncC, j = AttackECP(TrCode) 
if A == "Error":
    print("Error: {0}".format(B))
else:
    print("t-ECP = (A, B) for C(D', G) Dual:\n\nA(t + g) = A:")
    print(A)
    print("\nB(t + g) = B:")
    print(B)

###################### Проверка атаки
#j = 0
#P = D[j]
#D.remove(P)
#AA = codes.EvaluationAGCode(D, (t + g) * P)
#BB = codes.EvaluationAGCode(D, G - (t + g) * P)
#print()
#print(AA.generator_matrix())
#print()
#print(BB.generator_matrix())
#print(AA.generator_matrix() == A and BB.generator_matrix() == B)
