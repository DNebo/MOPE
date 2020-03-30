from random import *
from math import *
from _pydecimal import Decimal
import scipy.stats
p = 0.95
q = round(1 - p,2)  # Рівень значимості
m = 3
f1 = m - 1  # Ступені свободи
f2 = 4
f3 = f1 * f2
d = 4
N = 4
fisher = 0
x1min = -30
x1max = 20
x2min = -70
x2max = -10
x3min = -10
x3max = -40
xcpmax = (x1max + x2max + x3max) / 3
xcpmin = (x1min + x2min + x3min) / 3
ymax = int(200 + xcpmax)
ymin = int(200 + xcpmin)
ylist = []
ymed = []
znach = []
ydisp = []

print(ymax,ymin,)


def medium(n, l):
    k = "Yсереднє" + str(n) + " = ("
    for i in range(len(l) - 1):
        k += str(l[i]) + "+"
    k += str(len(l)) + ")/3 = " + str(ymed[n - 1])
    return k


def disp(n, l):
    global ydisp
    a = 0
    for i in range(len(l)):
        a += (l[i] - ymed[n - 1]) ** 2
    a /= len(l) - 1
    ydisp.append(round(a, 2))
    k = "D" + "{y" + str(n) + "} = ("
    for i in range(len(l) - 1):
        k += "(" + str(l[i]) + "-(" + str(ymed[n - 1]) + "))^2 + "
    k += "(" + str(l[len(l) - 1]) + "-(" + str(ymed[n - 1]) + "))^2)/" + str(len(l)) + " = " + str(ydisp[n - 1])
    return k


def equation():
    global ymed, ylist, a1, a2, a3, my, b

    for k in range(4):
        ylist.append([randrange(ymin, ymax) for i in range(m)])

    for k in range(4):
        ymed.append(round(sum(ylist[k]) / m, 2))
    my = round(sum(ymed) / 4, 2)  # Коефіціент b0 норм.

    a1 = round((-1 * (ymed[0] + ymed[1]) + 1 * (ymed[2] + ymed[3])) / 4, 2)  # Коефіціенти b1-3 норм.
    a2 = round((-1 * (ymed[0] + ymed[2]) + 1 * (ymed[1] + ymed[3])) / 4, 2)
    a3 = round((-1 * (ymed[0] + ymed[3]) + 1 * (ymed[1] + ymed[2])) / 4, 2)

    deltax1 = (fabs(x1max - x1min) / 2)
    deltax2 = (fabs(x2max - x2min) / 2)
    deltax3 = (fabs(x3max - x3min) / 2)

    x10 = (x1max + x1min) / 2
    x20 = (x2max + x2min) / 2
    x30 = (x3max + x3min) / 2

    b0 = round(my - (a1 * x10 / deltax1) - (a2 * x20 / deltax2) - (a3 * x30 / deltax3), 2)  # Коефіціент b0 натур.
    b1 = round(a1 / deltax1, 2)  # Коефіціенти b1-3 натур.
    b2 = round(a2 / deltax2, 2)
    b3 = round(a3 / deltax3, 2)
    b = [b0, b1, b2, b3]

    print("Генеруємо", m, "функцій відгуку для 3 експериметнів:\n Y1=", ylist[0], "\n Y2=", ylist[1], "\n Y3=",
          ylist[2],
          "\n Y3=", ylist[3], "\n")
    print("1)Знайдемо коефіціенти рівняння регресії:\n\n  Для нормованих значень:")
    print("  ", medium(1, ylist[0]), "\n  ", medium(2, ylist[1]), "\n  ", medium(3, ylist[2]), "\n  ",
          medium(4, ylist[3]))
    print("   b0 = my = ({0}+{1}+{2}+{3})/4 = {4}".format(*ymed, my))
    print("   b1 = a1 = (-1)*({0}+{1}) + 1*({2}+{3}) = {4}\n   b2 = a2 = (-1)*({0}+{2}) + 1*({1}+{3}) = {5}\n   \
b3 = a3 = (-1)*({0}+{3}) + 1*({1}+{2}) = {6}\n".format(*ymed, a1, a2, a3))
    print("  Для натуральних:")
    print("   Δx1 = |({0})-({1})|/2 = {2}".format(x1max, x1min, deltax1))
    print("   Δx2 = |({0})-({1})|/2 = {2}".format(x2max, x2min, deltax2))
    print("   Δx3 = |({0})-({1})|/2 = {2}".format(x3max, x3min, deltax3))
    print("   x10 = (({0})+({1}))/2 = {2}".format(x1max, x1min, x10))
    print("   x20 = (({0})+({1}))/2 = {2}".format(x2max, x2min, x20))
    print("   x30 = (({0})+({1}))/2 = {2}".format(x3max, x3min, x30))
    print("   b0 = {0} - ({1})*{2}/{3} - ({4})*{5}/{6} - ({4})*{5}/{6} - ({7})*{8}/{9} = {10}" \
          .format(my, a1, x10, deltax1, a2, x20, deltax2, a3, x30, deltax3, b0))
    print("   b1 = {0}\n   b2 = {1}\n   b3 = {2}".format(b1, b2, b3))


def Cochran():
    d1 = disp(1, ylist[0])
    d2 = disp(2, ylist[1])
    d3 = disp(3, ylist[2])
    d4 = disp(4, ylist[3])
    groz = round(max(ydisp) / sum(ydisp), 2)
    partresult = q / (f2 - 1)
    params = [partresult, f1, (f2 - 2) * f1]
    fisher = scipy.stats.f.isf(*params)
    result = fisher / (fisher + (f2 - 2))
    gkr = round(Decimal(result).quantize(Decimal('.0001')).__float__(), 2)

    print("\n2)Критерій Кохрана:\n\n  Знайдемо дисперсії по рядках:")
    print("  ", d1, "\n  ", d2, "\n  ", d3, "\n  ", d4, "\n")
    print("   Dmax{{yi}} = {0}\n   Gp = {0}/({1}+{2}+{3}+{4}) = {5}".format(max(ydisp), *ydisp, groz))
    print("   f1 = {0} - 1 = {1}, f2 = 4, q = {3}\n   За таблицею Gкр = {2}".format(m, f1, gkr, q))
    if groz < gkr:
        print("   Gp < Gкр => За критерієм Кохрана дисперсія однорідна з ймовірністю", p)
    else:
        print("   Gp > Gкр => За критерієм Кохрана дисперсія неоднорідна з ймовірністю", p)


def Student():
    global znach, y1, y2, y3, y4, disp, d
    disp = round(sum(ydisp) / 4, 2)
    dbs = round(d / (4 * m), 2)
    sbs = round(sqrt(dbs), 2)

    t0 = round(fabs(my) / sbs, 2)
    t1 = round(fabs(a1) / sbs, 2)
    t2 = round(fabs(a2) / sbs, 2)
    t3 = round(fabs(a3) / sbs, 2)
    tlist = [t0, t1, t2, t3]
    tkr = Decimal(abs(scipy.stats.t.ppf(q / 2, f3))).quantize(Decimal('.0001')).__float__()

    for troz in tlist:
        if troz < tkr:
            b[tlist.index(troz)] = 0
            d -= 1

    y1 = round(b[0] + b[1] * x1min + b[2] * x2min + b[3] * x3min, 2)
    y2 = round(b[0] + b[1] * x1min + b[2] * x2max + b[3] * x3max, 2)
    y3 = round(b[0] + b[1] * x1max + b[2] * x2min + b[3] * x3max, 2)
    y4 = round(b[0] + b[1] * x1max + b[2] * x2max + b[3] * x3min, 2)

    print("\n2)Критерій Стьюдента:\n")
    print("   Dвідтворюваності = ({0}+{1}+{2}+{3})/4 = {4}".format(*ydisp, disp))
    print("   D{{bi}} = {0}/(4*{1}) = {2}\n   S{{bi}} = sqrt({2}) = {3}".format(disp, m, dbs, sbs))
    print("   t0 = |{0}|/{4} = {5}\n   t1 = |{1}|/{4} = {6}\n   t2 = |{2}|/{4} = {7}\n   \
t3 = |{3}|/{4} = {8}\n   ".format(my, a1, a2, a3, sbs, *tlist))
    print("   f3 = 4*({0}-1) = {1}\n   За таблицею tkr = {2}".format(m, f3, tkr))
    for i in range(4):
        if tlist[i] < tkr:
            print(
                "   {0} < {1} => За критерієм Стьюдента коефіцієнт b{2} статистично незначущий з ймовірністю {3}".format(
                    tlist[i], tkr, i, p))
        else:
            print(
                "   {0} > {1} => За критерієм Стьюдента коефіцієнт b{2} статистично значимий з ймовірністю {3}".format(
                    tlist[i], tkr, i, p))
    print("\n   {0} + {1}*{4} + {2}*{5} + {3}*{6} = {7}".format(*b, x1min, x2min, x3min, y1))
    print("   {0} + {1}*{4} + {2}*{5} + {3}*{6} = {7}".format(*b, x1min, x2max, x3max, y2))
    print("   {0} + {1}*{4} + {2}*{5} + {3}*{6} = {7}".format(*b, x1max, x2min, x3max, y3))
    print("   {0} + {1}*{4} + {2}*{5} + {3}*{6} = {7}".format(*b, x1max, x2max, x3min, y4))


def Fisher():

      sad = round(m * ((y1 - ymed[0]) ** 2 + (y2 - ymed[1]) ** 2 + (y3 - ymed[2]) ** 2 + (y4 - ymed[3]) ** 2) / (4 - d),
                2)
      froz = round(sad / disp, 2)

      f4 = N - d

      fkr = Decimal(abs(scipy.stats.f.isf(q, f4, f3))).quantize(Decimal('.0001')).__float__()


      print("\n3)Критерій Фішера:\n")
      print("   f4 = {2} - {0} = {1}".format(d, f4, N))
      print(
          "   {0}*(({5} - {1})**2 + ({6} - {2})**2 + ({7} - {2})**2 + ({8} - {2})**2)/(4-{10}) = {9}".format(m, *ymed, y1,
                                                                                                           y2, y3,
                                                                                                           y4, sad, d))
      print("   Fр = {0}/{1} = {2}".format(sad, disp, froz))
      print("   За таблицею Fкр =", fkr)
      if fkr > froz:
        print("   За критерієм Фішера рівняння регресії адекватне оригіналу з ймовірністю", p)
      else:
        print("   За критерієм Фішера рівняння регресії неадекватне оригіналу з ймовірністю", p)


equation()
Cochran()
Student()
Fisher()

