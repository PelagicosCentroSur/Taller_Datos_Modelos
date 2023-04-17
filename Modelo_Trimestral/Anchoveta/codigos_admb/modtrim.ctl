#	CONTROLES_DEL_MODELO ANCHOVETA																												
#_________________________________________________________________________________																													
#	Coeficiente	de	variación	de	los	reclutamientos																							
#	sigmaR																												
0.6																													
#dt_desove  dt_indices       																										
0.5	    0.3         	     																										
#	PARAMETROS	DE	CRECIMIENTO	Y	MORTALIDAD	NATURAL																							
#____________________________________________________																													
#Loo	k	Lr	sr	b	M	h      																							
#20.8	0.11	6.0	1.0	0.5	0.22	0.75  #0.75   # M/k=2
20.8	0.11	2.36	1.0	0.5	0.22	1.0  #0.75   # M/k=2
0.1	0.3	10.1	0.05	0.2	0.1	0.1     #coeficientes	de variacion de	las prior (solo	activo	Lr)													

#	lag de tiempo con el reclutamiento (trimestres)																												
1																													
#____________________________________________________																													
#	Coeficiente_de_hiper_estabilidad/reducción_(si_es_1,_es_proporcional)																												
1																													
#	Parámetros_de_Selectividad																												
#____________________________________________________																													
#L50 (cm), s1	y s2 (modelo domo) (L50= talla	modal	app)																		
10	2	1000
#	opción_selectividad_(<0,logistico, >0, domo)_FLOTA																									
-1
#L50 pelas (cm), s_pel
10      1

#L50 reclas (cm), s_rec
10      1

#	Número_de_bloques_de_selectividad_(logistica)_y_años_de_inicio_FLOTA																												
1
1990.1


#	Número_de_bloques_de_capturabilidad_y_años_de_inicio_FLOTA																												
1
1990.1
#_________________________________________________________________																													
#	FASES_DE_ESTIMACION_DE_PARAMETROS																												

#	Fase_de_estimacion_de_capturabilidad_Flota																												
6       																												
#	Fase_estimación_de_la_selectividad_flota																												
3       																												
#	Fase_estimación_de_la_selectividad_y capturabilidad del PELACES 																												
4   5
#	Fase_estimación_de_la_selectividad_y capturabilidad del RECLAS																												
4   4
#	Fase_de_estimacion_de_Lr_sr_beta																											
5      -4	-6																											
#	Fase_de_estimacion_mortalidad_por_pesca																												
2																													
#	Fase_de_estimacion_variación_anual_Rt y Estacional_Rs																												
2	2																												
#	Fase_de_estimacion_desvios_condición_inicial_(si <0 esta en equilibrio)																								
-3																													
#	Fase_de_estimacion_parámetro_de_hiper_estabilidad_de_la_CPUE																												
-3																													

#	PROYECCION_DE_LA_POBLACION																												
#____________________________________________________																													
#numero_casos_de multiplicadores de F																													
13																											
#
1.00E-05	0.1	0.2	0.3	0.4	0.5	0.6	0.7	0.8	0.9	1	1.2	1.5																													
#	Años_a_proyectar																												
10																									
