# thermal_systems_project1

Problem Statement 

As a coolant for an automobile engine, Ethylene Glycol (50% diluted with deionized water) exits the engine at 90°C and flows through a heat exchanger that should cool the fluid to 60°C. The coolant flows at 1.2kg/s and is cooled using ambient air at 30°C. The Glycol side is susceptible to fouling, and the flow inlet and outlet headers of the heat exchanger have been ignored for the purpose of this design.


Optimization Techniques 

Initially, certain values were randomly changed, and patterns were observed. The first condition to be altered was the length of the tube. For the above conditions in section 4.2 to be achieved, the tube had to be 359.23m long, which is an unreasonable value for a car radiator.  

As such, the tube length was chosen to fit within range of 1m to 1.75 m, while other variables were changed. All calculations, optimization iterations and graphs were done using MATLAB R2022a. Within this fixed range for tube length, other dimensions on the heat exchanger were altered and the overall heat transfer coefficient, U, was derived. The parameters continuously changed including the longitudinal pitch (which directly affects the fin width), the transverse pitch (which directly affects the fin length), fin thickness, number of tubes, number of fins, number of rows and tube diameter. The physics derived U value was then compared to the expected U value. After each alteration, the percentage difference between the derived U and 
Uexpected was computed for the tube lengths of 1m to 1.75m and plotted with tube length on the x-axis and percentage difference between derived U and expected U on the y-axis. The point at which the graph cuts the x-axis shows the tube-length where U = Uexpected, thus satisfying our first condition. If the graph did not cross the x-axis, our first condition was not satisfied within tube lengths of 1m to 1.75m. This design failure can be observed in Figure 3 below. Hence, another set of alterations was carried out on our variable design parameters until U = Uexpected was satisfied. 


Qualitative Discussion of Pressure Drop and Pumping Power 

Since fluids are pumped through heat exchangers, it is important to know the amount of pumping power required for the system to perform at optimum rates as per design specifications. Another reason is that this determines the operational costs, and this advises the feasibility of the design proposal. Pressure drop is directly related to pumping power because as the pressure drop increases, the required power to pump the fluid increases as well.  

The proposed model has a low pumping power requirement of 631.932W. This suggests that our proposed model consumes less power to operate and therefore is cost-efficient.  

  