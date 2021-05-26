import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spint

matrix= np.array([[1,2,3],[4,5,6]])


def ManualGeometryDefinition():
    Tab_wezlow = np.array([[1, 0],
                      [2, 1],
                      [3, 0.5],
                      [4, 0.75]])

    Tab_Elementow = np.array([[1, 1, 3],
                      [2, 4, 2],
                      [3, 3, 4]])

    War_Brzeg = [{"ind": 1, "typ": 'D', "wartosc": 1},
          {"ind": 2, "typ": 'D', "wartosc": 2}]

    return Tab_wezlow, Tab_Elementow, War_Brzeg


def ShowGeometry (tab_w): # utworzona przeze mnie funkcja rusujaca geometrie
    y = np.zeros(tab_w.shape[0])
    plt.plot(tab_w[:, 1], y, marker='o') #wyrysowanie elementow
    for i in range(0, np.size(y), 1): #podpis wezlow
        plt.text(x=tab_w[i, 1], y=y[i]-0.007, s=int(tab_w[i, 0]),fontsize=7,color ='green' )
        # plt.text(x=tab_w[i, 1]+, y=y[i] - 0.007, s=tab_w[i, 0], fontsize=12, color='green')
    for i in range(0, np.size(y) - 1, 1): #podpis elementow
        plt.text(x=(tab_w[i, 1] + tab_w[i + 1, 1]) / 2, y=y[i] + 0.003, s=int(i + 1),fontsize=7, color='blue')

def AutomaticGeometryDefinition (x_a,x_b, n): # utworzona przeze mnie funkcja generujaca automatycznie geometrie
    # l_elementow = n-1
    odstep = (x_b - x_a) / (n - 1) #wyliczenie odstepu
    wezly = np.array([1, x_a]) #pierwszy wiersz tablicy wezly
    elementy = np.array([1,1,2]) #pierwszy wiersz tablicy elementy
    for i in range(1, n, 1):
       wezly = np.block([
            [wezly ],
            [i+1, i * odstep+x_a],
        ])
    for j in range(2, n,1):
        elementy = np.block([
            [elementy],
            [j,j,j+1]])
    return wezly,elementy



def Allocation(n):

    A = np.zeros([n, n]) #macierz glonbalna
    b = np.zeros([n, 1]) # wektor prawostronny

    return A, b


def BaseFunctions(stop_fun_baz):

    if stop_fun_baz == 0:
        funkcja = (lambda x: 1 + 0 * x)
        pochodna = (lambda x: 0 * x)

    elif stop_fun_baz == 1:

        funkcja = (lambda x: -1 / 2 * x + 1 / 2, lambda x: 0.5 * x + 0.5)
        pochodna = (lambda x: -1 / 2 + 0 * x, lambda x: 0.5 + 0 * x)

    elif stop_fun_baz == 2:
        funkcja =(lambda x: 0.5* x * (x - 1), lambda x: -x ** 2 + 1, lambda x: 0.5 * x * (x + 1))
        pochodna =(lambda x: x - 0.5, lambda x: -2 * x, lambda x: x + 0.5)

    else:
        raise Exception("Nieoczekiwany stop_fun_baz w BaseFunctions().")

    return funkcja, pochodna


def Aij(pochodna_fun_i, pochodna_fun_j, c, funkcja_i, funkcja_j):
    # funkcja podcałkowa
    funkcja_podcalkowa = lambda x: -pochodna_fun_i(x) * pochodna_fun_j(x) + c * funkcja_i(x) * funkcja_j(x)

    return funkcja_podcalkowa

def ShowSolution (Tab_wezlow, rozwiazanie):

    ShowGeometry(Tab_wezlow)

    x = Tab_wezlow[:, 1]

    plt.plot(x, rozwiazanie, 'm*')
    plt.show()


if __name__ == '__main__':

    # Preprocessing

    ## parametry sterujace
    c = 0
    f = lambda x: 0 * x  # wymuszenie

    ## 2. Geometria

    # zadefiniowana przez reczne utworzenie tabel

    # WEZLY, ELEMENTY, WB = ManualGeometryDefinition()
    # n = np.shape(WEZLY)[0]

    #Automatyczne wygenerowanie geometrii

    x_a = 0
    x_b = 1
    n = 5
    WEZLY, ELEMENTY = AutomaticGeometryDefinition(x_a, x_b, n)
    # warunki brzegowe
    WB = [{"ind": 1, "typ": 'D', "wartosc": 1},
          {"ind": n, "typ": 'D', "wartosc": 2}]


    #print(WEZLY)
    #print('\n')
    #print(ELEMENTY)
    #print('\n')
    #print(WB)

    #ShowGeometry(WEZLY)
    #plt.show()

    A, b = Allocation(n)

    # print(A)
    # print(b)

    stopien_fun_bazowych = 1      # 0-2
    phi, dphi = BaseFunctions(stopien_fun_bazowych)

    # x = np.linspace(-1,1, 101)
    # plt.plot(x, phi[0](x), 'red')
    # plt.plot(x, phi[1](x), 'green')
    # plt.plot(x, phi[2](x), 'brown')
    # plt.plot(x, dphi[0](x), 'blue')
    # plt.plot(x, dphi[1](x), 'pink')
    # plt.plot(x, dphi[2](x), 'black')
    # plt.show()

    # PROCESSING

    ElemsNumber = np.shape(ELEMENTY)[0]

    for l in np.arange(0, ElemsNumber):
        Glob_Ind_ele = ELEMENTY[l, 0]
        ind_Wezel_pocz = ELEMENTY[l, 1]  # indeks wezla poczatkowego elemntu ee
        ind_Wezel_kon = ELEMENTY[l, 2]  # indeks wezla koncowego elemntu ee
        ind_Glob_Wezlow = np.array([ind_Wezel_pocz, ind_Wezel_kon])

        Ml = np.zeros([stopien_fun_bazowych + 1, stopien_fun_bazowych + 1])


        x_a = WEZLY[ind_Wezel_pocz - 1, 1]  #indeksy pythonowe
        x_b = WEZLY[ind_Wezel_kon - 1, 1]

        J = (x_b - x_a) / 2 #Jakobian

        m = 0
        n = 0
        Ml[m, n] = J * spint.quad(Aij(dphi[m], dphi[n], c, phi[m], phi[n]), -1, 1)[0]

        m = 0
        n = 1
        Ml[m, n] = J * spint.quad(Aij(dphi[m], dphi[n], c, phi[m], phi[n]), -1, 1)[0]

        m = 1
        n = 0
        Ml[m, n] = J * spint.quad(Aij(dphi[m], dphi[n], c, phi[m], phi[n]), -1, 1)[0]

        m = 1
        n = 1
        Ml[m, n] = J * spint.quad(Aij(dphi[m], dphi[n], c, phi[m], phi[n]), -1, 1)[0]

        A[np.ix_(ind_Glob_Wezlow - 1, ind_Glob_Wezlow - 1)] = \
            A[np.ix_(ind_Glob_Wezlow - 1, ind_Glob_Wezlow - 1)] + Ml

        # print(Ml)
        # print('\n')
    print(A)
    print('\n')
    '''
    print(A)
    print('\n')
    Aprobna=A
    Aprobna_wyciecie = np.delete(Aprobna,[0,4],0)
    Aprobna_wyciecie1 = Aprobna_wyciecie
    Aprobna_wyciecie1 = np.delete(Aprobna_wyciecie1,[0,4],1)
    print(Aprobna_wyciecie1)
    #print(WB)
    '''
    # UWZGLEDNIENIE WARUNKOW BRZEGOWYCH
    '''
    if WB[0]['typ'] == 'D':
        ind_wezla = WB[0]['ind']
        wart_war_brzeg = WB[0]['wartosc']

        iwp = ind_wezla - 1

        WZMACNIACZ = 10 ** 14

        b[iwp] = A[iwp, iwp] * WZMACNIACZ * wart_war_brzeg
        A[iwp, iwp] = A[iwp, iwp] * WZMACNIACZ
    '''
   # print(A)
    if (WB[0]['typ'] == 'D') & (WB[1]['typ'] == 'D'):

        ind_wezla_1 = WB[0]['ind']   # indeks globalny wezla pocz
        wart_war_brzeg_1 = WB[0]['wartosc'] # war wezla pocz

        ind_wezla_2 = WB[1]['ind'] # indeks globalny wezla kon
        wart_war_brzeg_2 = WB[1]['wartosc'] # war wezla kon

        iwp1 = ind_wezla_1 - 1 # indeksy zmienione na indeksy pythona
        iwp2 = ind_wezla_2 - 1

        A=np.delete(A,[iwp1,iwp2],0) # usuwanie wierszy z macierzy A


        b=np.delete(b,[iwp1,iwp2],0) # usuwanie wierszy z wektora b

        n_b = np.shape(b)[0] # wielkosc wektora b

        for i in np.arange(0, n_b): # uwzglednienie wpływu wartosci w warunkach brzegowych na zmiane wartosci wektora b
            b[i] = b[i] - A[i,iwp1] * wart_war_brzeg_1
            b[i] = b[i] - A[i, iwp2] * wart_war_brzeg_2

        A = np.delete(A, [iwp1, iwp2], 1) # usuwanie kolumn z wartosci A

    '''
    if WB[1]['typ'] == 'D':
        ind_wezla = WB[1]['ind']
        wart_war_brzeg = WB[1]['wartosc']

        iwp = ind_wezla - 1

        WZMACNIACZ = 10 ** 14

        b[iwp] = A[iwp, iwp] * WZMACNIACZ * wart_war_brzeg
        A[iwp, iwp] = A[iwp, iwp] * WZMACNIACZ
    '''
    if WB[0]['typ'] == 'N':
        print('Nie zaimplementowano jeszcze. Zad.dom')

    if WB[1]['typ'] == 'N':
        print('Nie zaimplementowano jeszcze. Zad.dom')


    print(A)
    print('\n')
    print(b)
    print('\n')
    # Rozwiazanie ukl row lin
    u = np.linalg.solve(A, b)
    u = np.vstack(([wart_war_brzeg_1], u, [wart_war_brzeg_2])) # dodanie wartosci brzegowych funkcji do rozwiazania
    print(u)

    ShowSolution(WEZLY, u)