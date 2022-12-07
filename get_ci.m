function ci = get_ci(alpha)
C_m	= [1.61194131476772	1.49393452012555	1.37528538869492	1.25610819255259	1.13651320707464	1.01660633645202	0.896488770728602	0.776256693653611	0.656001072093091	0.535807580252877	0.415756761336965	0.295924647741176	0.176384363074465	0.057210030001782	-0.061513514137235	-0.179665611292373	-0.29702194928402	-0.413057842597917	-0.526238124587493	-0.628576560240532	-0.646406854082793	-0.853970209806321	-0.990078283758672	-1.11578355307185	-1.23933553203202	-1.36219603186014	-1.48490266989545	-1.60766791880069	-1.7304914630866	-1.85348245907399	-1.97652656114249	-2.09952917207188	-2.22239295497099	-2.34499972381981	-2.46734627244827	-2.58922594281404	-2.7105198313849	-2.83113361613618	-2.95086615867188	-3.06959165426447	-3.18725916452806];
poly = polyfit([-20:20], C_m, 7);

ci = poly(1) * alpha.^7 + poly(2) * alpha.^6 + poly(3) * alpha.^5 + poly(4) * alpha.^4 + poly(5) * alpha.^3 + poly(6) * alpha.^2 + poly(7) * alpha + poly(8);

end