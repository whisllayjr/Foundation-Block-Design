import tkinter as tk
from tkinter import ttk 
import math as m
from PIL import ImageTk, Image

fyd = 5000/1.15

def converte_fck(fck):
	# MPa -> kgf/cm^2
	return fck*10

def calc_lado_bloco(diam_est):
	return 2*diam_est

def calc_lado_A(diam_est):
	return 5*diam_est

def calc_lado_B(diam_est):
	return 2*diam_est

def calc_ap_eq(ap, bp):
	return round(m.sqrt(ap*bp), 2)

def calc_d_1est(diam_est):
	return 2*diam_est
'''
def calc_d_234est(qtd_est, e, ap):
	if qtd_est == "2est":
		return round(e/2 - ap/4)
	elif qtd_est == "3est":
		return  round(e*m.sqrt(3)/3 - 0.3*ap)
	elif qtd_est == "4est":
		return round(e*m.sqrt(2)/2 - ap*m.sqrt(2)/4)
'''

def calc_d_2est(e, ap):
	# alpha=45
	return round(e/2 - ap/4)

def calc_d_min34est(qtd_est, e, ap):
	if qtd_est == "3est":
		return round(0.58*(e - ap/2), 2)
	elif qtd_est == "4est":
		return round(0.71*(e - ap/2), 2)


def calc_d_max34est(qtd_est, e, ap):
	if qtd_est == "3est":
		return round(0.825*(e - ap/2), 2)
	elif qtd_est == "4est":
		return round(e - ap/2, 2)

def calc_d_adotado(d_min, d_max):
	d_min = m.ceil(d_min)
	d_max = m.floor(d_max)

	dlist = [35, 40, 45, 50, 55, 60, 65]

	for d in dlist:
		if d > d_min:
			return d
		else:
			continue


def calc_e(qtd_est, diam_est):
	if qtd_est == "2est":
		return 3*diam_est

	elif qtd_est == ("3est" and "4est"):
		return 2.5*diam_est


def calc_alpha(qtd_est, d, e, ap):
	ap_eq = ap
	if qtd_est == "2est":
		return round(m.degrees(m.atan(d/(e/2 - ap/4))), 3)

	elif qtd_est == "3est":
		return round(m.degrees(m.atan(d/(e*m.sqrt(3)/3 - 0.3*ap_eq))), 3)

	elif qtd_est == "4est":
		return round(m.degrees(m.atan(d/(e*m.sqrt(2)/2 - ap_eq*m.sqrt(2)/4))), 3)


def calc_Rs(qtd_est, carga, e, d, ap_eq):
	N = carga*10
	if qtd_est == "3est":
		return round(N/9*((e*m.sqrt(3) - 0.9*ap_eq)/d), 2)
	elif qtd_est == "4est":
		return round((N*m.sqrt(2)*(2*e - ap_eq))/(16*d), 2)

def calc_R_linha_s(Rs):
	return round(Rs*m.sqrt(3)/3, 2)

def calc_Rc(qtd_est, carga, alpha):
	N = carga*10
	alpha = m.radians(alpha)
	if qtd_est == "3est":
		return round(N/(3*m.sin(alpha)), 2)
	elif qtd_est == "4est":
		return round(N/(4*m.sin(alpha)), 2)

def calc_Re_max(Nk, My, e):
	Nk *= 10
	My *= 10
	e /= 100
	Re_max = round(1.02*Nk/4 + 0.5*My/e, 2)
	return Re_max/10

def verifica_Re_max(Re_max, capacidade_est):
	if Re_max <= capacidade_est:
		return "ok!"
	else:
		return "não ok!"

def calc_Nd(carga):
	# tf = 10*kN
	return carga*10*1.4

def calc_area_p(ap, bp):
	return ap*bp


def calc_sigma_cd_b_pil(qtd_est, Nd, area_p, alpha):
	if qtd_est == ("2est" and "3est"):
		return round(Nd / (area_p * (m.sin(m.radians(alpha))**2)), 2)
	elif qtd_est == "4est":
		return round(Nd / (area_p * (m.sin(m.radians(alpha))**2)), 2)

def calc_sigma_cd_b_est(qtd_est, Nd, diam_est, alpha):
	area_e = m.pi*diam_est**2/4
	if qtd_est == "2est":
		return round(Nd / (2*area_e * (m.sin(m.radians(alpha))**2)), 2)

	elif qtd_est == "3est":
		return round(Nd / (3*area_e * (m.sin(m.radians(alpha))**2)), 2)

	elif qtd_est == "4est":
		return round(Nd / (4*area_e * (m.sin(m.radians(alpha))**2)), 2)



def calc_tensao_lim(qtd_est, fck):
	# MPa -> kN/cm²
	fck /=10
	fcd = fck/1.4
	K_R = 0.9
	if qtd_est == "2est":
		return round(1.4*K_R*fcd, 2)

	elif qtd_est == "3est":
		return round(1.75*K_R*fcd, 2)

	elif qtd_est == "4est":
		return round(2.1*K_R*fcd, 2)

def verifica_tensaolim(tensao_atuante, tensao_lim):
	if tensao_atuante < tensao_lim:
		return 'ok!'
	else:
		return 'não passou na verificação das bielas'

def verifica_alpha(alpha):
	if 40 <= alpha <= 55:
		return 'ok!'
	else:
		return 'alpha não passou na verificação'

def calc_As_eh(carga, d, lado_bloco, fyd):
	# tf -> kgf
	carga *= 1000

	return round((0.25*carga*1.4*d/lado_bloco)/fyd, 2)

def calc_Ac(lado_bloco):
	return lado_bloco**2

def calc_taxa_armadura(carga, Ac, fck, fyd):
	return round((2*carga*1.4 - Ac*fck)/(fyd*Ac), 4)

def verifica_taxa_armadura(taxa_armadura):
	if taxa_armadura < 0.2/100:
		return f'adotado valor min: 0,2%'
	else:
		return 'ok!'

def calc_As_ew(taxa_armadura, lado_bloco):
	return taxa_armadura * (lado_bloco**2)

def calc_As_barra(diam_barra):
	# mm -> cm
	diam_barra /= 10
	return (m.pi*diam_barra**2)/4

def calc_quant_estribo(qtd_armacoes, As, As_barra):
	if qtd_armacoes == "1armacao":
		N = As/As_barra
		N = str(N)
		ptdec = N.index(".")
		if N[ptdec+1] == "0":
			return round(float(N))
		elif N[ptdec+1] != "0":
			return m.ceil(float(N))

	elif qtd_armacoes == "2armacoes":
		return m.ceil(As/(2*As_barra))



def calc_As_princ(Nd, d, fyd, e, ap):
	fyd /= 100
	return round((1.15*Nd*(2*e - ap))/(8*d*fyd), 2)

def calc_As_sup(qtd_est, As_princ):
	if qtd_est == ("2est" and "3est"): 
		return round(0.2*As_princ, 2)
	elif qtd_est == "4est":
		return round(0.2*4*As_princ, 2)

def calc_As_pele(qtd_est, As_lado):
	if qtd_est == "3est":
		As_tot = 3*As_lado
		return round(As_tot/8, 2)
	elif qtd_est == "4est":
		As_tot = 4*As_lado
		return round(As_tot/8, 2)

def calc_As_pele_sw(b):
	return round(0.075*b, 4)

def calc_As_lado3est(R_linha_s, fyd):
	R_linha_sd = R_linha_s*1.4
	fyd /= 100
	return round(R_linha_sd/fyd, 2)

def calc_As_lado4est(Nd, d, fyd, e, ap):
	fyd /= 100
	return round(Nd*(2*e - ap)/(16*d*fyd), 2)

def calc_As_malha(qtd_est, As_lado):
	if qtd_est == "3est":
		return round(As_lado/5, 2)
	elif qtd_est == "4est":
		return round(0.25*As_lado, 2)


def calc_As_susp_tot(qtd_est, Nd, fyd):
	fyd /= 100
	if qtd_est == "3est":
		return round(Nd/(4.5*fyd), 2)
	elif qtd_est == "4est":
		return round(Nd/(6*fyd), 2)

def calc_As_susp_face(qtd_est, As_susp_tot):
	if qtd_est == "3est":
		return round(As_susp_tot/3, 2)
	elif qtd_est == "4est":
		return round(As_susp_tot/4, 2)


def calc_espacamento(quant_estribo):
	return round(100/quant_estribo)

#def detalha_As_princ()


# - - - - - - - - - - - - - - - - - - - - - - - - - -
# window GUI
window = tk.Tk()
window.title("Dimensionamento de Blocos")
window['padx'] = 5
window['pady'] = 5

# - - - - - - - - - - - - - - - - - - - - - - - - - -
# Notebook
notebook = ttk.Notebook(window)
notebook.grid(row=0, column=0, sticky=tk.E + tk.W + tk.N + tk.S, padx=30, pady=4)

# widgets
tab1 = tk.Frame(notebook)
tab2 = tk.Frame(notebook)
tab3 = tk.Frame(notebook)
tab4 = tk.Frame(notebook)

notebook.add(tab1, text="Bloco de 1 estaca", compound=tk.TOP)
notebook.add(tab2, text="Bloco de 2 estacas", compound=tk.TOP)
notebook.add(tab3, text="Bloco de 3 estacas", compound=tk.TOP)
notebook.add(tab4, text="Bloco de 4 estacas", compound=tk.TOP)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# TAB 1: BLOCO DE 1 ESTACA

# - - - - - - - - - - - - - - - - - - - - - - - - - -
# LabelFrame Dados Iniciais TAB 1
ditab1 = ttk.LabelFrame(tab1, text="Dados Iniciais")
ditab1.grid(row=0, column=0, padx=30, pady=4)

entry_f_tab1 = ttk.Frame(ditab1)
entry_f_tab1.grid(row=0, column=0, sticky="N")

img_f_tab1 = ttk.Frame(ditab1)
img_f_tab1.grid(row=0, column=1)

# widgets
carga_l_tab1 = ttk.Label(entry_f_tab1, text="carga P (tf)")
carga_l_tab1.grid(row=0, column=0, sticky="W")
carga_e_tab1 = ttk.Entry(entry_f_tab1, width=10)
carga_e_tab1.grid(row=0, column=1)

fck_l_tab1 = ttk.Label(entry_f_tab1, text="fck (MPa)")
fck_l_tab1.grid(row=1, column=0, sticky="W")
fck_e_tab1 = ttk.Entry(entry_f_tab1, width=10)
fck_e_tab1.grid(row=1, column=1, sticky="W")

diam_est_l_tab1 = ttk.Label(entry_f_tab1, text="diâmetro da estaca (cm)")
diam_est_l_tab1.grid(row=2, column=0, sticky="W")
diam_est_e_tab1 = ttk.Entry(entry_f_tab1, width=10)
diam_est_e_tab1.grid(row=2, column=1, sticky="W")

adotar_Aeh_l_tab1 = ttk.Label(entry_f_tab1, text="adote estribo horizontal (mm)")
adotar_Aeh_l_tab1.grid(row=3, column=0)
adotar_Aeh_cbx_value_tab1 = tk.StringVar()
adotar_Aeh_cbx_tab1 = ttk.Combobox(entry_f_tab1, height=4, width=7)
adotar_Aeh_cbx_tab1['textvariable'] = adotar_Aeh_cbx_value_tab1
adotar_Aeh_cbx_tab1.grid(row=3, column=1)
adotar_Aeh_cbx_tab1['values'] = ("6.3", "5")

adotar_Aew_l_tab1 = ttk.Label(entry_f_tab1, text="adote estribo vertical (mm)")
adotar_Aew_l_tab1.grid(row=4, column=0, sticky="W")
adotar_Aew_cbx_value_tab1 = tk.StringVar()
adotar_Aew_cbx_tab1 = ttk.Combobox(entry_f_tab1, height=4, width=7)
adotar_Aew_cbx_tab1['textvariable'] = adotar_Aew_cbx_value_tab1
adotar_Aew_cbx_tab1.grid(row=4, column=1)
adotar_Aew_cbx_tab1['values'] = ("10", "6.3")

# imagem ilustrativa DADOS INICIAIS TAB 1

ditab1_img = ImageTk.PhotoImage(Image.open(r"di_bloco1est.PNG"))
ditab1_img_l = ttk.Label(img_f_tab1, image=ditab1_img)
ditab1_img_l.grid(row=0, column=0, rowspan=10)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# LabelFrame Resultados TAB 1
rtab1 = ttk.LabelFrame(tab1, text="Resultados")
rtab1.grid(row=0, column=1, padx=30, pady=4, sticky="N"+"W")

# widgets
Aeh_l_tab1 = ttk.Label(rtab1, text="Aeh (cm^2) =")
Aeh_l_tab1.grid(row=0, column=0, sticky="W")

quant_estribo_l_h_tab1 = ttk.Label(rtab1, text="N hor=")
quant_estribo_l_h_tab1.grid(row=1, column=0, sticky="W")

taxa_armadura_l_tab1 = ttk.Label(rtab1, text="taxa_armadura =")
taxa_armadura_l_tab1.grid(row=2, column=0, sticky="W")

verificacaorho_l_tab1 = ttk.Label(rtab1, text="verificação rho min =")
verificacaorho_l_tab1.grid(row=3, column=0, sticky="W")

As_ew_l_tab1 = ttk.Label(rtab1, text="Aew (cm^2) =")
As_ew_l_tab1.grid(row=4, column=0, sticky="W")

quant_estribo_l_v_tab1 = ttk.Label(rtab1, text="N vert=")
quant_estribo_l_v_tab1.grid(row=5, column=0, sticky="W")

def button_calcular1est():

	carga = float(carga_e_tab1.get())
	fck = float(fck_e_tab1.get())
	#cnom = float(cnom_e_tab1.get())
	diam_est = float(diam_est_e_tab1.get())
	#d_linha = float(d_linha_e_tab1.get())

	# calculo
	fck = converte_fck(fck)
	lado_bloco = calc_lado_bloco(diam_est)
	d = calc_d_1est(diam_est)
	As_eh = calc_As_eh(carga, d, lado_bloco, fyd)
	As_eh_l_r = ttk.Label(rtab1, text=As_eh)
	As_eh_l_r.grid(row=0, column=1)
	
	diam_barra = float(adotar_Aeh_cbx_value_tab1.get())

	As_barra = calc_As_barra(diam_barra)
	quant_estribo = calc_quant_estribo("2armacoes", As_eh, As_barra)
	quant_estribo_l_r = ttk.Label(rtab1, text=quant_estribo)
	quant_estribo_l_r.grid(row=1, column=1)

	Ac = calc_Ac(lado_bloco)
	taxa_armadura = calc_taxa_armadura(carga, Ac, fck, fyd)
	taxa_armadura_l_r = ttk.Label(rtab1, text=taxa_armadura)
	taxa_armadura_l_r.grid(row=2, column=1)

	verificacaorho = verifica_taxa_armadura(taxa_armadura)
	verificacaorho_l_r = ttk.Label(rtab1, text=verificacaorho)
	verificacaorho_l_r.grid(row=3, column=1)

	if verificacaorho == 'adotado valor min: 0,2%':
		taxa_armadura = 0.2/100
		As_ew = calc_As_ew(taxa_armadura, lado_bloco)
		As_ew_l_r = ttk.Label(rtab1, text=As_ew)
		As_ew_l_r.grid(row=4, column=1)

		diam_barra = float(adotar_Aew_cbx_value_tab1.get())
		As_barra = calc_As_barra(diam_barra)
		quant_estribo = calc_quant_estribo("2armacoes", As_ew, As_barra)
		quant_estribo_l_r = ttk.Label(rtab1, text=quant_estribo)
		quant_estribo_l_r.grid(row=5, column=1)


calcbuttontab1 = ttk.Button(entry_f_tab1, text="calcular!")
calcbuttontab1.grid(row=5, column=0)
calcbuttontab1['command'] = button_calcular1est


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# TAB 2: BLOCO DE 2 ESTACAS

# - - - - - - - - - - - - - - - - - - - - - - - - - -
# LabelFrame Dados Iniciais TAB 2
ditab2 = ttk.LabelFrame(tab2, text="Dados Iniciais")
ditab2.grid(row=0, column=0, padx=30, pady=4, sticky="N")

# widgets
carga_l_tab2 = ttk.Label(ditab2, text="carga (tf)")
carga_l_tab2.grid(row=0, column=0, sticky="W")
carga_e_tab2 = ttk.Entry(ditab2, width=10)
carga_e_tab2.grid(row=0, column=1, sticky="W")

fck_l_tab2 = ttk.Label(ditab2, text="fck (MPa)")
fck_l_tab2.grid(row=1, column=0, sticky="W")
fck_e_tab2 = ttk.Entry(ditab2, width=10)
fck_e_tab2.grid(row=1, column=1, sticky="W")


diam_est_l_tab2 = ttk.Label(ditab2, text="diam_est (cm)")
diam_est_l_tab2.grid(row=3, column=0, sticky="W")
diam_est_e_tab2 = ttk.Entry(ditab2, width=10)
diam_est_e_tab2.grid(row=3, column=1, sticky="W")

ap_l_tab2 = ttk.Label(ditab2, text="ap (cm)")
ap_l_tab2.grid(row=5, column=0, sticky="W")
ap_e_tab2 = ttk.Entry(ditab2, width=10)
ap_e_tab2.grid(row=5, column=1, sticky="W")

bp_l_tab2 = ttk.Label(ditab2, text="bp (cm)")
bp_l_tab2.grid(row=6, column=0, sticky="W")
bp_e_tab2 = ttk.Entry(ditab2, width=10)
bp_e_tab2.grid(row=6, column=1, sticky="W")

adotar_As_princ_l = ttk.Label(ditab2, text="adote barra do As principal (mm)")
adotar_As_princ_l.grid(row=7, column=0, sticky="W")
adotar_As_princ_cbx_value = tk.StringVar()
adotar_As_princ_cbx = ttk.Combobox(ditab2, height=4, width=7)
adotar_As_princ_cbx['textvariable'] = adotar_As_princ_cbx_value
adotar_As_princ_cbx.grid(row=7, column=1, sticky="W")
adotar_As_princ_cbx['values'] = ("10", "12.5", "16")

adotar_As_sup_l_tab2 = ttk.Label(ditab2, text="adote barra do As sup (mm)")
adotar_As_sup_l_tab2.grid(row=8, column=0, sticky="W")
adotar_As_sup_cbx_value_tab2 = tk.StringVar()
adotar_As_sup_cbx_tab2 = ttk.Combobox(ditab2, height=4, width=7)
adotar_As_sup_cbx_tab2['textvariable'] = adotar_As_sup_cbx_value_tab2
adotar_As_sup_cbx_tab2.grid(row=8, column=1, sticky="W")
adotar_As_sup_cbx_tab2['values'] = ("8", "10", "12.5", "16")

adotar_As_pele_l_tab2 = ttk.Label(ditab2, text="adote barra do As pele (mm)")
adotar_As_pele_l_tab2.grid(row=9, column=0, sticky="W")
adotar_As_pele_cbx_value_tab2 = tk.StringVar()
adotar_As_pele_cbx_tab2 = ttk.Combobox(ditab2, height=4, width=7)
adotar_As_pele_cbx_tab2['textvariable'] = adotar_As_pele_cbx_value_tab2
adotar_As_pele_cbx_tab2.grid(row=9, column=1, sticky="W")
adotar_As_pele_cbx_tab2['values'] = ("8", "10", "12.5", "16")

adotar_Asw_l = ttk.Label(ditab2, text="adote barra do Asw (mm)")
adotar_Asw_l.grid(row=10, column=0, sticky="W")
adotar_Asw_cbx_value = tk.StringVar()
adotar_Asw_cbx = ttk.Combobox(ditab2, height=4, width=7)
adotar_Asw_cbx['textvariable'] = adotar_Asw_cbx_value
adotar_Asw_cbx.grid(row=10, column=1, sticky="W")
adotar_Asw_cbx['values'] = ("8", "10", "12.5", "16")

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# LabelFrame Resultados TAB 2
rtab2 = ttk.LabelFrame(tab2, text="Resultados")
rtab2.grid(row=0, column=1, padx=30, pady=4)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# widgets TAB 2
d_l = ttk.Label(rtab2, text="d (cm) =")
d_l.grid(row=0, column=0, sticky="W")

alpha_l_tab2 = ttk.Label(rtab2, text="alpha (º) =")
alpha_l_tab2.grid(row=1, column=0, sticky="W")

verificacaoalpha_l_tab2 = ttk.Label(rtab2, text="verificação alpha =")
verificacaoalpha_l_tab2.grid(row=2, column=0, sticky="W")

tensao_lim_l_tab2 = ttk.Label(rtab2, text="tensão lim (kN/cm² =")
tensao_lim_l_tab2.grid(row=3, column=0, sticky="W")

sigma_cd_b_pil_l_tab2 = ttk.Label(rtab2, text="tensão pilar (kN/cm^2)=")
sigma_cd_b_pil_l_tab2.grid(row=4, column=0, sticky="W")

sigma_cd_b_est_l_tab2 = ttk.Label(rtab2, text="tensão estaca (kN/cm^2)=")
sigma_cd_b_est_l_tab2.grid(row=5, column=0, sticky="W")

verificacaotensaolimpil_l_tab2 = ttk.Label(rtab2, text="verificação tensão lim pilar =")
verificacaotensaolimpil_l_tab2.grid(row=6, column=0, sticky="W")

verificacaotensaolimest_l_tab2 = ttk.Label(rtab2, text="verificação tensão lim estaca =")
verificacaotensaolimest_l_tab2.grid(row=7, column=0, sticky="W")

As_princ_l = ttk.Label(rtab2, text="As principal =")
As_princ_l.grid(row=8, column=0, sticky="W")

quant_estribo_As_princ_l = ttk.Label(rtab2, text="N principal = ")
quant_estribo_As_princ_l.grid(row=9, column=0, sticky="W")

As_sup_l_tab2 = ttk.Label(rtab2, text="As sup = ")
As_sup_l_tab2.grid(row=10, column=0, sticky="W")

quant_estribo_As_sup_l_tab2 = ttk.Label(rtab2, text="N sup = ")
quant_estribo_As_sup_l_tab2.grid(row=11, column=0, sticky="W")

As_pele_l_tab3 = ttk.Label(rtab2, text="As pele = ")
As_pele_l_tab3.grid(row=12, column=0, sticky="W")

quant_estribo_As_pele_l_tab3 = ttk.Label(rtab2, text="N pele = ")
quant_estribo_As_pele_l_tab3.grid(row=13, column=0, sticky="W")

espacamento_As_pele_l_tab3 = ttk.Label(rtab2, text="S As pele = ")
espacamento_As_pele_l_tab3.grid(row=14, column=0, sticky="W")

Asw_l = ttk.Label(rtab2, text="Asw = ")
Asw_l.grid(row=15, column=0, sticky="W")

quant_estribo_Asw_l = ttk.Label(rtab2, text="N vert = ")
quant_estribo_Asw_l.grid(row=16, column=0, sticky="W")

espacamento_Asw_l = ttk.Label(rtab2, text="S Asw = ")
espacamento_Asw_l.grid(row=17, column=0, sticky="W")

def button_calcular2est():

	carga = float(carga_e_tab2.get())
	fck = float(fck_e_tab2.get())
	diam_est = float(diam_est_e_tab2.get())
	ap = float(ap_e_tab2.get())
	bp = float(bp_e_tab2.get())

	# cálculo
	Nd = calc_Nd(carga)
	area_p = calc_area_p(ap, bp)
	e = calc_e("2est", diam_est)
	b = calc_lado_B(diam_est)

	d = calc_d_2est(e, ap)
	d_l_r = ttk.Label(rtab2, text=d)
	d_l_r.grid(row=0, column=1)

	alpha = calc_alpha("2est", d, e, ap)
	alpha_l_r = ttk.Label(rtab2, text=alpha)
	alpha_l_r.grid(row=1, column=1)

	verificacaoalpha = verifica_alpha(alpha)
	verificacaoalpha_l_r = ttk.Label(rtab2, text=verificacaoalpha)
	verificacaoalpha_l_r.grid(row=2, column=1)

	tensao_lim = calc_tensao_lim("2est", fck)
	tensao_lim_l_r = ttk.Label(rtab2, text=tensao_lim)
	tensao_lim_l_r.grid(row=3, column=1)

	sigma_cd_b_pil = calc_sigma_cd_b_pil("2est", Nd, area_p, alpha)
	sigma_cd_b_pil_l_r = ttk.Label(rtab2, text=sigma_cd_b_pil)
	sigma_cd_b_pil_l_r.grid(row=4, column=1)

	sigma_cd_b_est = calc_sigma_cd_b_est("2est", Nd, diam_est, alpha)
	sigma_cd_b_est_l_r = ttk.Label(rtab2, text=sigma_cd_b_est)
	sigma_cd_b_est_l_r.grid(row=5, column=1)

	verificacaotensaolimpil = verifica_tensaolim(sigma_cd_b_pil, tensao_lim)
	verificacaotensaolimpil_l_r = ttk.Label(rtab2, text=verificacaotensaolimpil)
	verificacaotensaolimpil_l_r.grid(row=6, column=1)

	verificacaotensaolimest = verifica_tensaolim(sigma_cd_b_pil, tensao_lim)
	verificacaotensaolimest_l_r = ttk.Label(rtab2, text=verificacaotensaolimest)
	verificacaotensaolimest_l_r.grid(row=7, column=1)

	As_princ = calc_As_princ(Nd, d, fyd, e, ap)
	As_princ_l_r = ttk.Label(rtab2, text=As_princ)
	As_princ_l_r.grid(row=8, column=1)

	diam_barra_As_princ = float(adotar_As_princ_cbx_value.get())
	As_barra_As_princ = calc_As_barra(diam_barra_As_princ)
	quant_estribo_As_princ = calc_quant_estribo("1armacao", As_princ, As_barra_As_princ)
	quant_estribo_As_princ_l_r = ttk.Label(rtab2, text=quant_estribo_As_princ)
	quant_estribo_As_princ_l_r.grid(row=9, column=1)

	As_sup = calc_As_sup("2est", As_princ)
	As_sup_l_r = ttk.Label(rtab2, text=As_sup)
	As_sup_l_r.grid(row=10, column=1)

	diam_barra_As_sup = float(adotar_As_sup_cbx_value_tab2.get())
	As_barra_As_sup = calc_As_barra(diam_barra_As_sup)
	quant_estribo_As_sup = calc_quant_estribo("1armacao", As_sup, As_barra_As_sup)
	quant_estribo_As_sup_l_r = ttk.Label(rtab2, text=quant_estribo_As_sup)
	quant_estribo_As_sup_l_r.grid(row=11, column=1)

	As_pele = calc_As_pele_sw(b)
	As_pele_l_r = ttk.Label(rtab2, text=As_pele)
	As_pele_l_r.grid(row=12, column=1)

	diam_barra_As_pele = float(adotar_As_pele_cbx_value_tab2.get())
	As_barra_As_pele = calc_As_barra(diam_barra_As_pele)
	quant_estribo_As_pele = calc_quant_estribo("1armacao", As_pele, As_barra_As_pele)
	quant_estribo_As_pele_l_r = ttk.Label(rtab2, text=quant_estribo_As_pele)
	quant_estribo_As_pele_l_r.grid(row=13, column=1)

	espacamento_As_pele = calc_espacamento(quant_estribo_As_pele)
	espacamento_As_pele_l_r = ttk.Label(rtab2, text=espacamento_As_pele)
	espacamento_As_pele_l_r.grid(row=14, column=1)

	Asw = calc_As_pele_sw(b)
	Asw_l_r = ttk.Label(rtab2, text=Asw)
	Asw_l_r.grid(row=15, column=1)

	diam_barra_Asw = float(adotar_Asw_cbx_value.get())
	As_barra_Asw = calc_As_barra(diam_barra_Asw)
	quant_estribo_Asw = calc_quant_estribo("1armacao", Asw, As_barra_Asw)
	quant_estribo_Asw_l_r = ttk.Label(rtab2, text=quant_estribo_Asw)
	quant_estribo_Asw_l_r.grid(row=16, column=1)

	espacamento_Asw = calc_espacamento(quant_estribo_Asw)
	espacamento_Asw_l_r = ttk.Label(rtab2, text=espacamento_Asw)
	espacamento_Asw_l_r.grid(row=17, column=1)
		
calcbuttontab2 = ttk.Button(ditab2, text="calcular!")
calcbuttontab2.grid(row=20, column=0)
calcbuttontab2['command'] = button_calcular2est

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# TAB 3: BLOCO DE 3 ESTACAS

# - - - - - - - - - - - - - - - - - - - - - - - - - -
# LabelFrame Dados Iniciais TAB 3
ditab3 = ttk.LabelFrame(tab3, text="Dados Iniciais")
ditab3.grid(row=0, column=0, padx=30, pady=4, sticky="N")

# widgets

carga_l_tab3 = ttk.Label(ditab3, text="carga (tf)")
carga_l_tab3.grid(row=0, column=0, sticky="W")
carga_e_tab3 = ttk.Entry(ditab3, width=10)
carga_e_tab3.grid(row=0, column=1, sticky="W")

fck_l_tab3 = ttk.Label(ditab3, text="fck (MPa)")
fck_l_tab3.grid(row=1, column=0, sticky="W")
fck_e_tab3 = ttk.Entry(ditab3, width=10)
fck_e_tab3.grid(row=1, column=1, sticky="W")

diam_est_l_tab3 = ttk.Label(ditab3, text="diam_est (cm)")
diam_est_l_tab3.grid(row=3, column=0, sticky="W")
diam_est_e_tab3 = ttk.Entry(ditab3, width=10)
diam_est_e_tab3.grid(row=3, column=1, sticky="W")

ap_l_tab3 = ttk.Label(ditab3, text="ap (cm)")
ap_l_tab3.grid(row=5, column=0, sticky="W")
ap_e_tab3 = ttk.Entry(ditab3, width=10)
ap_e_tab3.grid(row=5, column=1, sticky="W")

bp_l_tab3 = ttk.Label(ditab3, text="bp (cm)")
bp_l_tab3.grid(row=6, column=0, sticky="W")
bp_e_tab3 = ttk.Entry(ditab3, width=10)
bp_e_tab3.grid(row=6, column=1, sticky="W")

adotar_As_lado_l = ttk.Label(ditab3, text="adote barra do As lado (mm)")
adotar_As_lado_l.grid(row=7, column=0, sticky="W")
adotar_As_lado_cbx_value = tk.StringVar()
adotar_As_lado_cbx = ttk.Combobox(ditab3, height=4, width=7)
adotar_As_lado_cbx['textvariable'] = adotar_As_lado_cbx_value
adotar_As_lado_cbx.grid(row=7, column=1, sticky="W")
adotar_As_lado_cbx['values'] = ("10", "12.5", "16")


adotar_As_malha_l = ttk.Label(ditab3, text="adote barra do As malha (mm)")
adotar_As_malha_l.grid(row=8, column=0, sticky="W")
adotar_As_malha_cbx_value = tk.StringVar()
adotar_As_malha_cbx = ttk.Combobox(ditab3, height=4, width=7)
adotar_As_malha_cbx['textvariable'] = adotar_As_malha_cbx_value
adotar_As_malha_cbx.grid(row=8, column=1, sticky="W")
adotar_As_malha_cbx['values'] = ("6.3" "8", "10")

adotar_As_sup_l_tab3 = ttk.Label(ditab3, text="adote barra do As sup (mm)")
adotar_As_sup_l_tab3.grid(row=9, column=0, sticky="W")
adotar_As_sup_cbx_value_tab3 = tk.StringVar()
adotar_As_sup_cbx_tab3 = ttk.Combobox(ditab3, height=4, width=7)
adotar_As_sup_cbx_tab3['textvariable'] = adotar_As_sup_cbx_value_tab3
adotar_As_sup_cbx_tab3.grid(row=9, column=1, sticky="W")
adotar_As_sup_cbx_tab3['values'] = ("8", "10", "12.5", "16")

adotar_As_pele_l_tab3 = ttk.Label(ditab3, text="adote barra do As pele (mm)")
adotar_As_pele_l_tab3.grid(row=10, column=0, sticky="W")
adotar_As_pele_cbx_value_tab3 = tk.StringVar()
adotar_As_pele_cbx_tab3 = ttk.Combobox(ditab3, height=4, width=7)
adotar_As_pele_cbx_tab3['textvariable'] = adotar_As_pele_cbx_value_tab3
adotar_As_pele_cbx_tab3.grid(row=10, column=1, sticky="W")
adotar_As_pele_cbx_tab3['values'] = ("8", "10", "12.5", "16")

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# LabelFrame Resultados TAB 3
rtab3 = ttk.LabelFrame(tab3, text="Resultados")
rtab3.grid(row=0, column=1, padx=30, pady=4)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# widgets Resultados TAB 3

d_min_l = ttk.Label(rtab3, text="d_min (cm) =")
d_min_l.grid(row=0, column=0, sticky="W")

d_max_l = ttk.Label(rtab3, text="d_max (cm) =")
d_max_l.grid(row=1, column=0, sticky="W")

d_adotado_l = ttk.Label(rtab3, text="d adotado (cm) =")
d_adotado_l.grid(row=2, column=0, sticky="W")

alpha_l_tab3 = ttk.Label(rtab3, text="alpha (º) =")
alpha_l_tab3.grid(row=3, column=0, sticky="W")

verificacaoalpha_l_tab3 = ttk.Label(rtab3, text="verificação alpha =")
verificacaoalpha_l_tab3.grid(row=4, column=0, sticky="W")

Rs_l = ttk.Label(rtab3, text="Rs =")
Rs_l.grid(row=5, column=0, sticky="W")

Rc_l = ttk.Label(rtab3, text="Rc =")
Rc_l.grid(row=6, column=0, sticky="W")

tensao_lim_l_tab3 = ttk.Label(rtab3, text="tensão lim (kN/cm² =")
tensao_lim_l_tab3.grid(row=7, column=0, sticky="W")

sigma_cd_b_pil_l_tab3 = ttk.Label(rtab3, text="tensão pilar (kN/cm^2)=")
sigma_cd_b_pil_l_tab3.grid(row=8, column=0, sticky="W")

sigma_cd_b_est_l_tab3 = ttk.Label(rtab3, text="tensão estaca (kN/cm^2)=")
sigma_cd_b_est_l_tab3.grid(row=9, column=0, sticky="W")

verificacaotensaolimpil_l_tab3 = ttk.Label(rtab3, text="verificação tensão lim pilar =")
verificacaotensaolimpil_l_tab3.grid(row=10, column=0, sticky="W")

verificacaotensaolimest_l_tab3 = ttk.Label(rtab3, text="verificação tensão lim estaca =")
verificacaotensaolimest_l_tab3.grid(row=11, column=0, sticky="W")

R_linha_s_l = ttk.Label(rtab3, text="R_linha_s =")
R_linha_s_l.grid(row=12, column=0, sticky="W")

As_lado_l = ttk.Label(rtab3, text="As_lado =")
As_lado_l.grid(row=13, column=0, sticky="W")

quant_estribo_As_lado_l = ttk.Label(rtab3, text="N lado = ")
quant_estribo_As_lado_l.grid(row=14, column=0, sticky="W")

As_malha_l = ttk.Label(rtab3, text="As_malha =")
As_malha_l.grid(row=15, column=0, sticky="W")

As_susp_tot_l = ttk.Label(rtab3, text="As_susp_tot =")
As_susp_tot_l.grid(row=16, column=0, sticky="W")

As_susp_face_l = ttk.Label(rtab3, text="As_susp_face =")
As_susp_face_l.grid(row=17, column=0, sticky="W")

quant_estribo_As_malha_l = ttk.Label(rtab3, text="N malha = ")
quant_estribo_As_malha_l.grid(row=19, column=0, sticky="W")


As_sup_l_tab3 = ttk.Label(rtab3, text="As sup = ")
As_sup_l_tab3.grid(row=20, column=0, sticky="W")

quant_estribo_As_sup_l_tab3 = ttk.Label(rtab3, text="N sup = ")
quant_estribo_As_sup_l_tab3.grid(row=21, column=0, sticky="W")

As_pele_l_tab3 = ttk.Label(rtab3, text="As pele = ")
As_pele_l_tab3.grid(row=22, column=0, sticky="W")

quant_estribo_As_pele_l_tab3 = ttk.Label(rtab3, text="N pele = ")
quant_estribo_As_pele_l_tab3.grid(row=23, column=0, sticky="W")

def button_calcular3est():
	carga = float(carga_e_tab3.get())
	fck = float(fck_e_tab3.get())
	diam_est = float(diam_est_e_tab3.get())
	ap = float(ap_e_tab3.get())
	bp = float(bp_e_tab3.get())

	# cálculo
	Nd = calc_Nd(carga)
	area_p = calc_area_p(ap, bp)
	'''
	b = calc_lado_B(diam_est)
	'''
	e = calc_e("3est", diam_est)
	ap_eq = calc_ap_eq(ap, bp)
	
	d_min = calc_d_min34est("3est", e, ap_eq)
	d_min_l_r = ttk.Label(rtab3, text=d_min)
	d_min_l_r.grid(row=0, column=1)

	d_max = calc_d_max34est("3est", e, ap_eq)
	d_max_l_r = ttk.Label(rtab3, text=d_max)
	d_max_l_r.grid(row=1, column=1)

	d_adotado = calc_d_adotado(d_min, d_max)
	d_adotado_l_r = ttk.Label(rtab3, text=d_adotado)
	d_adotado_l_r.grid(row=2, column=1)

	alpha = calc_alpha("3est", d_adotado, e, ap_eq)
	alpha_l_r = ttk.Label(rtab3, text=alpha)
	alpha_l_r.grid(row=3, column=1)

	verificacaoalpha = verifica_alpha(alpha)
	verificacaoalpha_l_r = ttk.Label(rtab3, text=verificacaoalpha)
	verificacaoalpha_l_r.grid(row=4, column=1)

	Rs = calc_Rs("3est", carga, e, d_adotado, ap_eq)
	Rs_l_r = ttk.Label(rtab3, text=Rs)
	Rs_l_r.grid(row=5, column=1)

	Rc = calc_Rc("3est", carga,alpha)
	Rc_l_r = ttk.Label(rtab3, text=Rc)
	Rc_l_r.grid(row=6, column=1)

	tensao_lim = calc_tensao_lim("3est", fck)
	tensao_lim_l_r = ttk.Label(rtab3, text=tensao_lim)
	tensao_lim_l_r.grid(row=7, column=1)

	sigma_cd_b_pil = calc_sigma_cd_b_pil("3est", Nd, area_p, alpha)
	sigma_cd_b_pil_l_r = ttk.Label(rtab3, text=sigma_cd_b_pil)
	sigma_cd_b_pil_l_r.grid(row=8, column=1)

	sigma_cd_b_est = calc_sigma_cd_b_est("3est", Nd, diam_est, alpha)
	sigma_cd_b_est_l_r = ttk.Label(rtab3, text=sigma_cd_b_est)
	sigma_cd_b_est_l_r.grid(row=9, column=1)

	verificacaotensaolimpil = verifica_tensaolim(sigma_cd_b_pil, tensao_lim)
	verificacaotensaolimpil_l_r = ttk.Label(rtab3, text=verificacaotensaolimpil)
	verificacaotensaolimpil_l_r.grid(row=10, column=1)

	verificacaotensaolimest = verifica_tensaolim(sigma_cd_b_pil, tensao_lim)
	verificacaotensaolimest_l_r = ttk.Label(rtab3, text=verificacaotensaolimest)
	verificacaotensaolimest_l_r.grid(row=11, column=1)

	R_linha_s = calc_R_linha_s(Rs)
	R_linha_s_l_r = ttk.Label(rtab3, text=R_linha_s)
	R_linha_s_l_r.grid(row=12, column=1)

	As_lado = calc_As_lado3est(R_linha_s, fyd)
	As_lado_l_r = ttk.Label(rtab3, text=As_lado)
	As_lado_l_r.grid(row=13, column=1)

	diam_barra_As_lado = float(adotar_As_lado_cbx_value.get())
	As_barra_As_lado = calc_As_barra(diam_barra_As_lado)
	quant_estribo_As_lado = calc_quant_estribo("1armacao", As_lado, As_barra_As_lado)
	quant_estribo_As_lado_l_r = ttk.Label(rtab3, text=quant_estribo_As_lado)
	quant_estribo_As_lado_l_r.grid(row=14, column=1)

	As_malha = calc_As_malha("3est", As_lado)
	As_malha_l_r = ttk.Label(rtab3, text=As_malha)
	As_malha_l_r.grid(row=15, column=1)

	As_susp_tot = calc_As_susp_tot("3est", Nd, fyd)
	As_susp_tot_l_r = ttk.Label(rtab3, text=As_susp_tot)
	As_susp_tot_l_r.grid(row=16, column=1)

	As_susp_face = calc_As_susp_face("3est", As_susp_tot)
	As_susp_face_l_r = ttk.Label(rtab3, text=As_susp_face)
	As_susp_face_l_r.grid(row=17, column=1)

	if As_susp_face >= As_malha:
		As_malha_l = ttk.Label(rtab3, text="As_malha = As_susp_face =")
		As_malha_l.grid(row=18, column=0)
		As_malha_l_r = ttk.Label(rtab3, text=As_susp_face)
		As_malha_l_r.grid(row=18, column=1)

		diam_barra_As_malha = float(adotar_As_malha_cbx_value.get())
		As_barra_As_malha = calc_As_barra(diam_barra_As_malha)
		quant_estribo_As_malha = calc_quant_estribo("1armacao", As_susp_face, As_barra_As_malha)
		quant_estribo_As_malha_l_r = ttk.Label(rtab3, text=quant_estribo_As_malha)
		quant_estribo_As_malha_l_r.grid(row=19, column=1)

	As_sup = calc_As_sup("3est", As_lado)
	As_sup_l_r = ttk.Label(rtab3, text=As_sup)
	As_sup_l_r.grid(row=20, column=1)

	diam_barra_As_sup = float(adotar_As_sup_cbx_value_tab3.get())
	As_barra_As_sup = calc_As_barra(diam_barra_As_sup)
	quant_estribo_As_sup = calc_quant_estribo("1armacao", As_sup, As_barra_As_sup)
	quant_estribo_As_sup_l_r = ttk.Label(rtab3, text=quant_estribo_As_sup)
	quant_estribo_As_sup_l_r.grid(row=21, column=1)

	As_pele = calc_As_pele("3est", As_lado)
	As_pele_l_r = ttk.Label(rtab3, text=As_pele)
	As_pele_l_r.grid(row=22, column=1)

	diam_barra_As_pele = float(adotar_As_pele_cbx_value_tab3.get())
	As_barra_As_pele = calc_As_barra(diam_barra_As_pele)
	quant_estribo_As_pele = calc_quant_estribo("1armacao", As_pele, As_barra_As_pele)
	quant_estribo_As_pele_l_r = ttk.Label(rtab3, text=quant_estribo_As_pele)
	quant_estribo_As_pele_l_r.grid(row=23, column=1)

calcbuttontab3 = ttk.Button(ditab3, text="calcular!")
calcbuttontab3.grid(row=20, column=0)
calcbuttontab3['command'] = button_calcular3est


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# TAB 4: BLOCO DE 4 ESTACAS

# - - - - - - - - - - - - - - - - - - - - - - - - - -
# LabelFrame Dados Iniciais TAB 4
ditab4 = ttk.LabelFrame(tab4, text="Dados Iniciais")
ditab4.grid(row=0, column=0, padx=30, pady=4, sticky="N")

# widgets Dados Iniciais TAB 4

carga_l_tab4 = ttk.Label(ditab4, text="carga (tf)")
carga_l_tab4.grid(row=0, column=0, sticky="W")
carga_e_tab4 = ttk.Entry(ditab4, width=10)
carga_e_tab4.grid(row=0, column=1, sticky="W")

My_l_tab4 = ttk.Label(ditab4, text="My (tf.m)")
My_l_tab4.grid(row=1, column=0, sticky="W")
My_e_tab4 = ttk.Entry(ditab4, width=10)
My_e_tab4.grid(row=1, column=1, sticky="W")

fck_l_tab4 = ttk.Label(ditab4, text="fck (MPa)")
fck_l_tab4.grid(row=2, column=0, sticky="W")
fck_e_tab4 = ttk.Entry(ditab4, width=10)
fck_e_tab4.grid(row=2, column=1, sticky="W")

diam_est_l_tab4 = ttk.Label(ditab4, text="diam_est (cm)")
diam_est_l_tab4.grid(row=3, column=0, sticky="W")
diam_est_e_tab4 = ttk.Entry(ditab4, width=10)
diam_est_e_tab4.grid(row=3, column=1, sticky="W")

capacidade_est_l_tab4 = ttk.Label(ditab4, text="capacidade_est (tf)")
capacidade_est_l_tab4.grid(row=4, column=0, sticky="W")
capacidade_est_e_tab4 = ttk.Entry(ditab4, width=10)
capacidade_est_e_tab4.grid(row=4, column=1, sticky="W")

adotar_carga_l = ttk.Label(ditab4, text="adote tipo de carregamento")
adotar_carga_l.grid(row=5, column=0, sticky="W")
adotar_carga_cbx_value = tk.StringVar()
adotar_carga_cbx = ttk.Combobox(ditab4, height=4, width=29)
adotar_carga_cbx['textvariable'] = adotar_carga_cbx_value
adotar_carga_cbx.grid(row=5, column=1, sticky="W")
adotar_carga_cbx['values'] = ("com carga de projeto: 4*Re_max", "com carga máxima das estacas")

ap_l_tab4 = ttk.Label(ditab4, text="ap (cm)")
ap_l_tab4.grid(row=6, column=0, sticky="W")
ap_e_tab4 = ttk.Entry(ditab4, width=10)
ap_e_tab4.grid(row=6, column=1, sticky="W")

bp_l_tab4 = ttk.Label(ditab4, text="bp (cm)")
bp_l_tab4.grid(row=7, column=0, sticky="W")
bp_e_tab4 = ttk.Entry(ditab4, width=10)
bp_e_tab4.grid(row=7, column=1, sticky="W")

adotar_As_lado_l = ttk.Label(ditab4, text="adote barra do As lado (mm)")
adotar_As_lado_l.grid(row=8, column=0, sticky="W")
adotar_As_lado_cbx_value = tk.StringVar()
adotar_As_lado_cbx = ttk.Combobox(ditab4, height=4, width=7)
adotar_As_lado_cbx['textvariable'] = adotar_As_lado_cbx_value
adotar_As_lado_cbx.grid(row=8, column=1, sticky="W")
adotar_As_lado_cbx['values'] = ("10", "12.5", "16")

adotar_As_malha_l = ttk.Label(ditab4, text="adote barra do As malha (mm)")
adotar_As_malha_l.grid(row=9, column=0, sticky="W")
adotar_As_malha_cbx_value = tk.StringVar()
adotar_As_malha_cbx = ttk.Combobox(ditab4, height=4, width=7)
adotar_As_malha_cbx['textvariable'] = adotar_As_malha_cbx_value
adotar_As_malha_cbx.grid(row=9, column=1, sticky="W")
adotar_As_malha_cbx['values'] = ("6.3", "8", "10")

adotar_As_sup_l_tab4 = ttk.Label(ditab4, text="adote barra do As sup (mm)")
adotar_As_sup_l_tab4.grid(row=10, column=0, sticky="W")
adotar_As_sup_cbx_value_tab4 = tk.StringVar()
adotar_As_sup_cbx_tab4 = ttk.Combobox(ditab4, height=4, width=7)
adotar_As_sup_cbx_tab4['textvariable'] = adotar_As_sup_cbx_value_tab4
adotar_As_sup_cbx_tab4.grid(row=10, column=1, sticky="W")
adotar_As_sup_cbx_tab4['values'] = ("8", "10", "12.5", "16")

adotar_As_pele_l_tab4 = ttk.Label(ditab4, text="adote barra do As pele (mm)")
adotar_As_pele_l_tab4.grid(row=11, column=0, sticky="W")
adotar_As_pele_cbx_value_tab4 = tk.StringVar()
adotar_As_pele_cbx_tab4 = ttk.Combobox(ditab4, height=4, width=7)
adotar_As_pele_cbx_tab4['textvariable'] = adotar_As_pele_cbx_value_tab4
adotar_As_pele_cbx_tab4.grid(row=11, column=1, sticky="W")
adotar_As_pele_cbx_tab4['values'] = ("8", "10", "12.5", "16")

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# LabelFrame Resultados TAB 4
rtab4 = ttk.LabelFrame(tab4, text="Resultados")
rtab4.grid(row=0, column=1, padx=30, pady=4)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# widgets Resultados TAB 4

d_min_l = ttk.Label(rtab4, text="d_min (cm) =")
d_min_l.grid(row=0, column=0, sticky="w")

d_max_l = ttk.Label(rtab4, text="d_max (cm) =")
d_max_l.grid(row=1, column=0, sticky="w")

d_adotado_l = ttk.Label(rtab4, text="d adotado (cm) =")
d_adotado_l.grid(row=2, column=0, sticky="w")

alpha_l_tab4 = ttk.Label(rtab4, text="alpha (º) =")
alpha_l_tab4.grid(row=3, column=0, sticky="w")

verificacaoalpha_l_tab4 = ttk.Label(rtab4, text="verificação alpha =")
verificacaoalpha_l_tab4.grid(row=4, column=0, sticky="w")

Re_max_l = ttk.Label(rtab4, text="Re_max =")
Re_max_l.grid(row=5, column=0, sticky="w")

verificacaoRe_max_l_tab4 = ttk.Label(rtab4, text="verificação Re_max =")
verificacaoRe_max_l_tab4.grid(row=6, column=0, sticky="w")

Rs_l = ttk.Label(rtab4, text="Rs =")
Rs_l.grid(row=7, column=0, sticky="w")

Rc_l = ttk.Label(rtab4, text="Rc =")
Rc_l.grid(row=8, column=0, sticky="w")

tensao_lim_l_tab4 = ttk.Label(rtab4, text="tensão lim (kN/cm² =")
tensao_lim_l_tab4.grid(row=9, column=0, sticky="w")

sigma_cd_b_pil_l_tab4 = ttk.Label(rtab4, text="tensão pilar (kN/cm^2)=")
sigma_cd_b_pil_l_tab4.grid(row=10, column=0, sticky="w")

sigma_cd_b_est_l_tab4 = ttk.Label(rtab4, text="tensão estaca (kN/cm^2)=")
sigma_cd_b_est_l_tab4.grid(row=11, column=0, sticky="w")

verificacaotensaolimpil_l_tab4 = ttk.Label(rtab4, text="verificação tensão lim pilar =")
verificacaotensaolimpil_l_tab4.grid(row=12, column=0, sticky="w")

verificacaotensaolimest_l_tab4 = ttk.Label(rtab4, text="verificação tensão lim estaca =")
verificacaotensaolimest_l_tab4.grid(row=13, column=0, sticky="w")

As_lado_l = ttk.Label(rtab4, text="As_lado =")
As_lado_l.grid(row=14, column=0, sticky="w")

quant_estribo_As_lado_l = ttk.Label(rtab4, text="N lado = ")
quant_estribo_As_lado_l.grid(row=15, column=0, sticky="w")

As_malha_l = ttk.Label(rtab4, text="As_malha =")
As_malha_l.grid(row=16, column=0, sticky="w")

As_susp_tot_l = ttk.Label(rtab4, text="As_susp_tot =")
As_susp_tot_l.grid(row=17, column=0, sticky="w")

As_susp_face_l = ttk.Label(rtab4, text="As_susp_face =")
As_susp_face_l.grid(row=18, column=0, sticky="w")

quant_estribo_As_malha_l = ttk.Label(rtab4, text="N malha = ")
quant_estribo_As_malha_l.grid(row=20, column=0, sticky="w")


As_sup_l_tab4 = ttk.Label(rtab4, text="As sup = ")
As_sup_l_tab4.grid(row=21, column=0, sticky="w")

quant_estribo_As_sup_l_tab4 = ttk.Label(rtab4, text="N sup = ")
quant_estribo_As_sup_l_tab4.grid(row=22, column=0, sticky="w")

As_pele_l_tab4 = ttk.Label(rtab4, text="As pele = ")
As_pele_l_tab4.grid(row=23, column=0, sticky="w")

quant_estribo_As_pele_l_tab4 = ttk.Label(rtab4, text="N pele = ")
quant_estribo_As_pele_l_tab4.grid(row=24, column=0, sticky="w")

def button_calcular4est():
	carga = float(carga_e_tab4.get())
	My = float(My_e_tab4.get())
	fck = float(fck_e_tab4.get())
	diam_est = float(diam_est_e_tab4.get())
	capacidade_est = float(capacidade_est_e_tab4.get())
	ap = float(ap_e_tab4.get())
	bp = float(bp_e_tab4.get())

	# cálculo
	Nd = calc_Nd(carga)
	area_p = calc_area_p(ap, bp)
	
	#b = calc_lado_B(diam_est)
	
	e = calc_e("4est", diam_est)
	ap_eq = calc_ap_eq(ap, bp)

	d_min = calc_d_min34est("4est", e, ap_eq)
	d_min_l_r = ttk.Label(rtab4, text=d_min)
	d_min_l_r.grid(row=0, column=1)

	d_max = calc_d_max34est("4est", e, ap_eq)
	d_max_l_r = ttk.Label(rtab4, text=d_max)
	d_max_l_r.grid(row=1, column=1)

	d_adotado = calc_d_adotado(d_min, d_max)
	d_adotado_l_r = ttk.Label(rtab4, text=d_adotado)
	d_adotado_l_r.grid(row=2, column=1)

	alpha = calc_alpha("4est", d_adotado, e, ap_eq)
	alpha_l_r = ttk.Label(rtab4, text=alpha)
	alpha_l_r.grid(row=3, column=1)

	verificacaoalpha = verifica_alpha(alpha)
	verificacaoalpha_l_r = ttk.Label(rtab4, text=verificacaoalpha)
	verificacaoalpha_l_r.grid(row=4, column=1)

	Re_max = calc_Re_max(carga, My, e)
	Re_max_l_r = ttk.Label(rtab4, text=Re_max)
	Re_max_l_r.grid(row=5, column=1)

	verificacaoRe_max = verifica_Re_max(Re_max, capacidade_est)
	verificacaoRe_max_l_r = ttk.Label(rtab4, text=verificacaoRe_max)
	verificacaoRe_max_l_r.grid(row=6, column=1)

	tipo_carregamento = adotar_carga_cbx_value.get()

	if tipo_carregamento == "com carga de projeto: 4*Re_max":
		carga = 4*Re_max
		Nd = calc_Nd(carga)

		Rs = calc_Rs("4est", carga, e, d_adotado, ap_eq)
		Rs_l_r = ttk.Label(rtab4, text=Rs)
		Rs_l_r.grid(row=7, column=1)

		Rc = calc_Rc("4est", carga, alpha)
		Rc_l_r = ttk.Label(rtab4, text=Rc)
		Rc_l_r.grid(row=8, column=1)

		tensao_lim = calc_tensao_lim("4est", fck)
		tensao_lim_l_r = ttk.Label(rtab4, text=tensao_lim)
		tensao_lim_l_r.grid(row=9, column=1)

		sigma_cd_b_pil = calc_sigma_cd_b_pil("4est", Nd, area_p, alpha)
		sigma_cd_b_pil_l_r = ttk.Label(rtab4, text=sigma_cd_b_pil)
		sigma_cd_b_pil_l_r.grid(row=10, column=1)

		sigma_cd_b_est = calc_sigma_cd_b_est("4est", Nd, diam_est, alpha)
		sigma_cd_b_est_l_r = ttk.Label(rtab4, text=sigma_cd_b_est)
		sigma_cd_b_est_l_r.grid(row=11, column=1)

		verificacaotensaolimpil = verifica_tensaolim(sigma_cd_b_pil, tensao_lim)
		verificacaotensaolimpil_l_r = ttk.Label(rtab4, text=verificacaotensaolimpil)
		verificacaotensaolimpil_l_r.grid(row=12, column=1)

		verificacaotensaolimest = verifica_tensaolim(sigma_cd_b_pil, tensao_lim)
		verificacaotensaolimest_l_r = ttk.Label(rtab4, text=verificacaotensaolimest)
		verificacaotensaolimest_l_r.grid(row=13, column=1)

		As_lado = calc_As_lado4est(Nd, d_adotado, fyd, e, ap_eq)
		As_lado_l_r = ttk.Label(rtab4, text=As_lado)
		As_lado_l_r.grid(row=14, column=1)

		diam_barra_As_lado = float(adotar_As_lado_cbx_value.get())
		As_barra_As_lado = calc_As_barra(diam_barra_As_lado)
		quant_estribo_As_lado = calc_quant_estribo("1armacao", As_lado, As_barra_As_lado)
		quant_estribo_As_lado_l_r = ttk.Label(rtab4, text=quant_estribo_As_lado)
		quant_estribo_As_lado_l_r.grid(row=15, column=1)

		As_malha = calc_As_malha("4est", As_lado)
		As_malha_l_r = ttk.Label(rtab4, text=As_malha)
		As_malha_l_r.grid(row=16, column=1)

		As_susp_tot = calc_As_susp_tot("4est", Nd, fyd)
		As_susp_tot_l_r = ttk.Label(rtab4, text=As_susp_tot)
		As_susp_tot_l_r.grid(row=17, column=1)

		As_susp_face = calc_As_susp_face("4est", As_susp_tot)
		As_susp_face_l_r = ttk.Label(rtab4, text=As_susp_face)
		As_susp_face_l_r.grid(row=18, column=1)

		if As_susp_face >= As_malha:
			As_malha_l = ttk.Label(rtab4, text="As_malha = As_susp_face =")
			As_malha_l.grid(row=19, column=0, sticky="W")
			As_malha_l_r = ttk.Label(rtab4, text=As_susp_face)
			As_malha_l_r.grid(row=19, column=1)

			diam_barra_As_malha = float(adotar_As_malha_cbx_value.get())
			As_barra_As_malha = calc_As_barra(diam_barra_As_malha)
			quant_estribo_As_malha = calc_quant_estribo("1armacao", As_susp_face, As_barra_As_malha)
			quant_estribo_As_malha_l_r = ttk.Label(rtab4, text=quant_estribo_As_malha)
			quant_estribo_As_malha_l_r.grid(row=20, column=1)

		As_sup = calc_As_sup("4est", As_lado)
		As_sup_l_r = ttk.Label(rtab4, text=As_sup)
		As_sup_l_r.grid(row=21, column=1)

		As_sup /= 2

		diam_barra_As_sup = float(adotar_As_sup_cbx_value_tab4.get())
		As_barra_As_sup = calc_As_barra(diam_barra_As_sup)
		quant_estribo_As_sup = calc_quant_estribo("1armacao", As_sup, As_barra_As_sup)
		quant_estribo_As_sup_l_r = ttk.Label(rtab4, text=quant_estribo_As_sup)
		quant_estribo_As_sup_l_r.grid(row=22, column=1)

		As_pele = calc_As_pele("4est", As_lado)
		As_pele_l_r = ttk.Label(rtab4, text=As_pele)
		As_pele_l_r.grid(row=23, column=1)

		diam_barra_As_pele = float(adotar_As_pele_cbx_value_tab4.get())
		As_barra_As_pele = calc_As_barra(diam_barra_As_pele)
		quant_estribo_As_pele = calc_quant_estribo("1armacao", As_pele, As_barra_As_pele)
		quant_estribo_As_pele_l_r = ttk.Label(rtab4, text=quant_estribo_As_pele)
		quant_estribo_As_pele_l_r.grid(row=24, column=1)

calcbuttontab4 = ttk.Button(ditab4, text="calcular!")
calcbuttontab4.grid(row=20, column=0)
calcbuttontab4['command'] = button_calcular4est

window.mainloop()
