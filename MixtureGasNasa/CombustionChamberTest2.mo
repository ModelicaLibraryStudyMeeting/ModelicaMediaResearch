package CombustionChamberTest2
  import SI = Modelica.SIunits;
  import ThermoPower.Choices;

  partial model CombustionChamberBase "Combustion Chamber"
    extends ThermoPower.Icons.Gas.Mixer;
    replaceable package Air = Modelica.Media.Interfaces.PartialMedium;
    replaceable package Fuel = Modelica.Media.Interfaces.PartialMedium;
    replaceable package Exhaust = Modelica.Media.Interfaces.PartialMedium;
    parameter SI.Volume V "Inner volume";
    parameter SI.Area S = 0 "Inner surface";
    parameter SI.CoefficientOfHeatTransfer gamma = 0 "Heat Transfer Coefficient" annotation(
      Evaluate = true);
    parameter SI.HeatCapacity Cm = 0 "Metal Heat Capacity" annotation(
      Evaluate = true);
    parameter SI.Temperature Tmstart = 300 "Metal wall start temperature" annotation(
      Dialog(tab = "Initialisation"));
    parameter SI.SpecificEnthalpy HH "Lower Heating value of fuel";
    parameter Boolean allowFlowReversal = system.allowFlowReversal "= true to allow flow reversal, false restricts to design direction" annotation(
      Evaluate = true);
    outer ThermoPower.System system "System wide properties";
    parameter Air.AbsolutePressure pstart = 101325 "Pressure start value" annotation(
      Dialog(tab = "Initialisation"));
    parameter Air.Temperature Tstart = 300 "Temperature start value" annotation(
      Dialog(tab = "Initialisation"));
    parameter Air.MassFraction Xstart[Exhaust.nX] = Exhaust.reference_X "Start flue gas composition" annotation(
      Dialog(tab = "Initialisation"));
    parameter Choices.Init.Options initOpt = system.initOpt "Initialisation option" annotation(
      Dialog(tab = "Initialisation"));
    parameter Boolean noInitialPressure = false "Remove initial equation on pressure" annotation(
      Dialog(tab = "Initialisation"),
      choices(checkBox = true));
    Exhaust.BaseProperties fluegas(p(start = pstart), T(start = Tstart), Xi(start = Xstart[1:Exhaust.nXi]));
    SI.Mass M "Gas total mass";
    SI.Mass MX[Exhaust.nXi] "Partial flue gas masses";
    SI.InternalEnergy E "Gas total energy";
    SI.Temperature Tm(start = Tmstart) "Wall temperature";
    Air.SpecificEnthalpy hia "Air specific enthalpy";
    Fuel.SpecificEnthalpy hif "Fuel specific enthalpy";
    Exhaust.SpecificEnthalpy ho "Outlet specific enthalpy";
    SI.Power HR "Heat rate";
    SI.Time Tr "Residence time";
    ThermoPower.Gas.FlangeA ina(redeclare package Medium = Air, m_flow(min = if allowFlowReversal then -Modelica.Constants.inf else 0)) "inlet air" annotation(
      Placement(transformation(extent = {{-120, -20}, {-80, 20}}, rotation = 0)));
    ThermoPower.Gas.FlangeA inf(redeclare package Medium = Fuel, m_flow(min = if allowFlowReversal then -Modelica.Constants.inf else 0)) "inlet fuel" annotation(
      Placement(transformation(extent = {{-20, 80}, {20, 120}}, rotation = 0)));
    ThermoPower.Gas.FlangeB out(redeclare package Medium = Exhaust, m_flow(max = if allowFlowReversal then +Modelica.Constants.inf else 0)) "flue gas" annotation(
      Placement(transformation(extent = {{80, -20}, {120, 20}}, rotation = 0)));
  equation
    M = fluegas.d * V "Gas mass";
    E = fluegas.u * M "Gas energy";
    MX = fluegas.Xi * M "Component masses";
    HR = inf.m_flow * HH;
    der(M) = ina.m_flow + inf.m_flow + out.m_flow "Gas mass balance";
    der(E) = ina.m_flow * hia + inf.m_flow * hif + out.m_flow * ho + HR - gamma * S * (fluegas.T - Tm) "Gas energy balance";
    if Cm > 0 and gamma > 0 then
      Cm * der(Tm) = gamma * S * (fluegas.T - Tm) "Metal wall energy balance";
    else
      Tm = fluegas.T;
    end if;
// Set gas properties
    out.p = fluegas.p;
    out.h_outflow = fluegas.h;
    out.Xi_outflow = fluegas.Xi;
// Boundary conditions
    ina.p = fluegas.p;
    ina.h_outflow = 0;
    ina.Xi_outflow = Air.reference_X[1:Air.nXi];
    inf.p = fluegas.p;
    inf.h_outflow = 0;
    inf.Xi_outflow = Fuel.reference_X[1:Fuel.nXi];
    assert(ina.m_flow >= 0, "The model does not support flow reversal");
    hia = inStream(ina.h_outflow);
    assert(inf.m_flow >= 0, "The model does not support flow reversal");
    hif = inStream(inf.h_outflow);
    assert(out.m_flow <= 0, "The model does not support flow reversal");
    ho = fluegas.h;
    Tr = noEvent(M / max(abs(out.m_flow), Modelica.Constants.eps));
  initial equation
// Initial conditions
    if initOpt == Choices.Init.Options.noInit then
// do nothing
    elseif initOpt == Choices.Init.Options.fixedState then
      if not noInitialPressure then
        fluegas.p = pstart;
      end if;
      fluegas.T = Tstart;
      fluegas.Xi = Xstart[1:Exhaust.nXi];
    elseif initOpt == Choices.Init.Options.steadyState then
      if not noInitialPressure then
        der(fluegas.p) = 0;
      end if;
      der(fluegas.T) = 0;
      der(fluegas.Xi) = zeros(Exhaust.nXi);
      if Cm > 0 and gamma > 0 then
        der(Tm) = 0;
      end if;
    elseif initOpt == Choices.Init.Options.steadyStateNoP then
      der(fluegas.T) = 0;
      der(fluegas.Xi) = zeros(Exhaust.nXi);
      if Cm > 0 and gamma > 0 then
        der(Tm) = 0;
      end if;
    else
      assert(false, "Unsupported initialisation option");
    end if;
    annotation(
      Documentation(info = "<html>
This is the model-base of a Combustion Chamber, with a constant volume.
<p>The metal wall temperature and the heat transfer coefficient between the wall and the fluid are uniform. The wall is thermally insulated from the outside. It has been assumed that inlet gases are premixed before entering in the volume.
<p><b>Modelling options</b></p>
<p>This model has three different Medium models to characterize the inlet air, fuel, and flue gas exhaust.
<p>If <tt>gamma = 0</tt>, the thermal effects of the surrounding walls are neglected.</p>
<p>There are two ways to obtain correct energy balances. The first is to explicitly set the lower heating value of the fuel <tt>HH</tt>, and use medium models that do not include the enthalpy of formation, by setting <tt>excludeEnthalpyOfFormation = true</tt>, which is the default option in Modelica.Media. As the heating value is usually provided at 25 degC temperature, it is also necessary to set <tt>referenceChoice =ReferenceEnthalpy.ZeroAt25C</tt> in all medium models for consistency. This is done in the medium models contained within <a href=\"modelica://ThermoPower.Media\">ThermoPower.Media</a>.</p>
<p>Alternatively, one can set <tt>excludeEnthalpyOfFormation = false</tt> in all media and set <tt>HH = 0</tt>. By doing so, the heating value is automatically accounted for by the difference in the enthalpy of formation. 
</html>", revisions = "<html>
<ul>
<li><i>30 May 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Initialisation support added.</li>
<li><i>31 Jan 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
    CombustionChamber model restructured using inheritance.
<p> First release.</li>
</ul>
</html>
      "),
      Diagram(graphics));
  end CombustionChamberBase;

  model CombustionChamber "Combustion Chamber"
    extends CombustionChamberBase(redeclare package Air = ThermoPower.Media.Air "O2, H2O, Ar, N2", redeclare package Fuel = ThermoPower.Media.NaturalGas "N2, CO2, CH4", redeclare package Exhaust = ThermoPower.Media.FlueGas "O2, Ar, H2O, CO2, N2");
    Real wcomb(final quantity = "MolarFlowRate", unit = "mol/s") "Molar Combustion rate (CH4)";
    SI.PerUnit lambda "Stoichiometric ratio (>1 if air flow is greater than stoichiometric)";
  protected
    Air.MassFraction ina_X[Air.nXi] = inStream(ina.Xi_outflow);
    Fuel.MassFraction inf_X[Fuel.nXi] = inStream(inf.Xi_outflow);
  equation
    wcomb = inf.m_flow * inf_X[3] / Fuel.data[3].MM "Combustion molar flow rate";
    lambda = ina.m_flow * ina_X[1] / Air.data[1].MM / (2 * wcomb);
    assert(lambda >= 1, "Not enough oxygen flow");
    der(MX[1]) = ina.m_flow * ina_X[1] + out.m_flow * fluegas.X[1] - 2 * wcomb * Exhaust.data[1].MM "oxygen";
    der(MX[2]) = ina.m_flow * ina_X[3] + out.m_flow * fluegas.X[2] "argon";
    der(MX[3]) = ina.m_flow * ina_X[2] + out.m_flow * fluegas.X[3] + 2 * wcomb * Exhaust.data[3].MM "water";
    der(MX[4]) = inf.m_flow * inf_X[2] + out.m_flow * fluegas.X[4] + wcomb * Exhaust.data[4].MM "carbondioxide";
    der(MX[5]) = ina.m_flow * ina_X[4] + out.m_flow * fluegas.X[5] + inf.m_flow * inf_X[1] "nitrogen";
    annotation(
      Icon(graphics),
      Documentation(info = "<html>
This model extends the CombustionChamber Base model, with the definition of the gases.
<p>In particular, the air inlet uses the <tt>Media.Air</tt> medium model, the fuel input uses the <tt>Media.NaturalGas</tt> medium model, and the flue gas outlet uses the <tt>Medium.FlueGas</tt> medium model.
<p>The composition of the outgoing gas is determined by the mass balance of every component, taking into account the combustion reaction CH4+2O2--->2H2O+CO2.</p>
<p>The model assumes complete combustion, so that it is only valid if the oxygen flow at the air inlet is greater than the stoichiometric flow corresponding to the flow at the fuel inlet.</p>

</html>", revisions = "<html>
<ul>
<li><i>31 Jan 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
 Combustion Chamber model restructured using inheritance.
     <p>  First release.
 </li>
</ul>
</html>"));
  end CombustionChamber;

  model TestCC1
    extends Modelica.Icons.Example;
    ThermoPower.Gas.SourceMassFlow Wcompressor(redeclare package Medium = ThermoPower.Media.Air, w0 = 158, T = 616.95) annotation(
      Placement(transformation(extent = {{-80, -10}, {-60, 10}}, rotation = 0)));
    CombustionChamber CombustionChamber1(initOpt = ThermoPower.Choices.Init.Options.steadyState, HH = 41.6e6, pstart = 11.2e5, V = 0.1, S = 0.1) annotation(
      Placement(transformation(extent = {{-38, -10}, {-18, 10}}, rotation = 0)));
    ThermoPower.Gas.SourceMassFlow Wfuel(redeclare package Medium = ThermoPower.Media.NaturalGas, use_in_w0 = true) annotation(
      Placement(transformation(extent = {{-50, 28}, {-30, 48}}, rotation = 0)));
    ThermoPower.Gas.PressDrop PressDrop1(redeclare package Medium = ThermoPower.Media.FlueGas, FFtype = ThermoPower.Choices.PressDrop.FFtypes.OpPoint, rhonom = 3.3, wnom = 158.9, pstart = 11.2e5, dpnom = 0.426e5) annotation(
      Placement(transformation(extent = {{-4, -10}, {16, 10}}, rotation = 0)));
    ThermoPower.Gas.SensT SensT1(redeclare package Medium = ThermoPower.Media.FlueGas) annotation(
      Placement(transformation(extent = {{26, -6}, {46, 14}}, rotation = 0)));
    Modelica.Blocks.Sources.Step Step1(startTime = 0.5, height = -0.3, offset = 3.1) annotation(
      Placement(transformation(extent = {{-78, 56}, {-58, 76}}, rotation = 0)));
    ThermoPower.Gas.ValveLin ValveLin1(redeclare package Medium = ThermoPower.Media.FlueGas, Kv = 161.1 / 9.77e5) annotation(
      Placement(transformation(extent = {{54, -10}, {74, 10}}, rotation = 0)));
    ThermoPower.Gas.SinkPressure SinkP1(redeclare package Medium = ThermoPower.Media.FlueGas) annotation(
      Placement(transformation(extent = {{84, -10}, {104, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Constant Constant1(k = 1) annotation(
      Placement(transformation(extent = {{22, 28}, {42, 48}}, rotation = 0)));
    inner ThermoPower.System system annotation(
      Placement(transformation(extent = {{80, 80}, {100, 100}})));
  equation
    connect(Wfuel.flange, CombustionChamber1.inf) annotation(
      Line(points = {{-30, 38}, {-28, 38}, {-28, 10}}, color = {159, 159, 223}, thickness = 0.5));
    connect(Wcompressor.flange, CombustionChamber1.ina) annotation(
      Line(points = {{-60, 0}, {-38, 0}}, color = {159, 159, 223}, thickness = 0.5));
    connect(CombustionChamber1.out, PressDrop1.inlet) annotation(
      Line(points = {{-18, 0}, {-4, 0}}, color = {159, 159, 223}, thickness = 0.5));
    connect(PressDrop1.outlet, SensT1.inlet) annotation(
      Line(points = {{16, 0}, {30, 0}}, color = {159, 159, 223}, thickness = 0.5));
    connect(Step1.y, Wfuel.in_w0) annotation(
      Line(points = {{-57, 66}, {-46, 66}, {-46, 43}}, color = {0, 0, 127}));
    connect(ValveLin1.outlet, SinkP1.flange) annotation(
      Line(points = {{74, 0}, {84, 0}}, color = {159, 159, 223}, thickness = 0.5));
    connect(SensT1.outlet, ValveLin1.inlet) annotation(
      Line(points = {{42, 0}, {54, 0}}, color = {159, 159, 223}, thickness = 0.5));
    connect(Constant1.y, ValveLin1.cmd) annotation(
      Line(points = {{43, 38}, {64, 38}, {64, 7}}, color = {0, 0, 127}));
    annotation(
      Documentation(info = "<html>
This model tests the <tt>CombustionChamber</tt> model. The model start at steady state. At time t = 0.5, the fuel flow rate is reduced by 10%.

<p>Simulate for 5s.
</html>"),
      experiment(StopTime = 5));
  end TestCC1;

  model TestCC2
    extends Modelica.Icons.Example;
    ThermoPower.Gas.SourceMassFlow Wcompressor(redeclare package Medium = ThermoPower.Media.Air, w0 = 158, T = 616.95) annotation(
      Placement(transformation(extent = {{-80, -10}, {-60, 10}}, rotation = 0)));
    CombustionChamber CombustionChamber1(HH = 0, S = 0.1, V = 0.1, initOpt = ThermoPower.Choices.Init.Options.steadyState, pstart = 11.2e5) annotation(
      Placement(transformation(extent = {{-38, -10}, {-18, 10}}, rotation = 0)));
    ThermoPower.Gas.SourceMassFlow Wfuel(redeclare package Medium = ThermoPower.Media.NaturalGas, use_in_w0 = true) annotation(
      Placement(transformation(extent = {{-50, 28}, {-30, 48}}, rotation = 0)));
    ThermoPower.Gas.PressDrop PressDrop1(redeclare package Medium = ThermoPower.Media.FlueGas, FFtype = ThermoPower.Choices.PressDrop.FFtypes.OpPoint, rhonom = 3.3, wnom = 158.9, pstart = 11.2e5, dpnom = 0.426e5) annotation(
      Placement(transformation(extent = {{-4, -10}, {16, 10}}, rotation = 0)));
    ThermoPower.Gas.SensT SensT1(redeclare package Medium = ThermoPower.Media.FlueGas) annotation(
      Placement(transformation(extent = {{26, -6}, {46, 14}}, rotation = 0)));
    Modelica.Blocks.Sources.Step Step1(startTime = 0.5, height = -0.3, offset = 3.1) annotation(
      Placement(transformation(extent = {{-78, 56}, {-58, 76}}, rotation = 0)));
    ThermoPower.Gas.ValveLin ValveLin1(redeclare package Medium = ThermoPower.Media.FlueGas, Kv = 161.1 / 9.77e5) annotation(
      Placement(transformation(extent = {{54, -10}, {74, 10}}, rotation = 0)));
    ThermoPower.Gas.SinkPressure SinkP1(redeclare package Medium = ThermoPower.Media.FlueGas) annotation(
      Placement(transformation(extent = {{84, -10}, {104, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Constant Constant1(k = 1) annotation(
      Placement(transformation(extent = {{22, 28}, {42, 48}}, rotation = 0)));
    inner ThermoPower.System system annotation(
      Placement(transformation(extent = {{80, 80}, {100, 100}})));
  equation
    connect(Wfuel.flange, CombustionChamber1.inf) annotation(
      Line(points = {{-30, 38}, {-28, 38}, {-28, 10}}, color = {159, 159, 223}, thickness = 0.5));
    connect(Wcompressor.flange, CombustionChamber1.ina) annotation(
      Line(points = {{-60, 0}, {-38, 0}}, color = {159, 159, 223}, thickness = 0.5));
    connect(CombustionChamber1.out, PressDrop1.inlet) annotation(
      Line(points = {{-18, 0}, {-4, 0}}, color = {159, 159, 223}, thickness = 0.5));
    connect(PressDrop1.outlet, SensT1.inlet) annotation(
      Line(points = {{16, 0}, {30, 0}}, color = {159, 159, 223}, thickness = 0.5));
    connect(Step1.y, Wfuel.in_w0) annotation(
      Line(points = {{-57, 66}, {-46, 66}, {-46, 43}}, color = {0, 0, 127}));
    connect(ValveLin1.outlet, SinkP1.flange) annotation(
      Line(points = {{74, 0}, {84, 0}}, color = {159, 159, 223}, thickness = 0.5));
    connect(SensT1.outlet, ValveLin1.inlet) annotation(
      Line(points = {{42, 0}, {54, 0}}, color = {159, 159, 223}, thickness = 0.5));
    connect(Constant1.y, ValveLin1.cmd) annotation(
      Line(points = {{43, 38}, {64, 38}, {64, 7}}, color = {0, 0, 127}));
    annotation(
      Documentation(info = "<html>
This model tests the <tt>CombustionChamber</tt> model. The model start at steady state. At time t = 0.5, the fuel flow rate is reduced by 10%.

<p>Simulate for 5s.
</html>"),
      experiment(StopTime = 5));
  end TestCC2;

  package Air "Air as mixture of O2, N2, Ar and H2O"
    extends Modelica.Media.IdealGases.Common.MixtureGasNasa(mediumName = "Air", data = {Modelica.Media.IdealGases.Common.SingleGasesData.O2, Modelica.Media.IdealGases.Common.SingleGasesData.H2O, Modelica.Media.IdealGases.Common.SingleGasesData.Ar, Modelica.Media.IdealGases.Common.SingleGasesData.N2}, fluidConstants = {Modelica.Media.IdealGases.Common.FluidData.O2, Modelica.Media.IdealGases.Common.FluidData.H2O, Modelica.Media.IdealGases.Common.FluidData.Ar, Modelica.Media.IdealGases.Common.FluidData.N2}, substanceNames = {"Oxygen", "Water", "Argon", "Nitrogen"}, reference_X = {0.23, 0.015, 0.005, 0.75}, excludeEnthalpyOfFormation = false, referenceChoice = Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C);
  end Air;

  package NaturalGas "Mixture of N2, CO2, and CH4"
    extends Modelica.Media.IdealGases.Common.MixtureGasNasa(mediumName = "NaturalGas", data = {Modelica.Media.IdealGases.Common.SingleGasesData.N2, Modelica.Media.IdealGases.Common.SingleGasesData.CO2, Modelica.Media.IdealGases.Common.SingleGasesData.CH4}, fluidConstants = {Modelica.Media.IdealGases.Common.FluidData.N2, Modelica.Media.IdealGases.Common.FluidData.CO2, Modelica.Media.IdealGases.Common.FluidData.CH4}, substanceNames = {"Nitrogen", "Carbondioxide", "Methane"}, reference_X = {0.02, 0.012, 0.968}, referenceChoice = Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C, excludeEnthalpyOfFormation = false);
  end NaturalGas;

  model TestCC3
    extends Modelica.Icons.Example;
    ThermoPower.Gas.SourceMassFlow Wcompressor(redeclare package Medium = Air, w0 = 158, T = 616.95) annotation(
      Placement(transformation(extent = {{-80, -10}, {-60, 10}}, rotation = 0)));
    CombustionChamber2 CombustionChamber1(HH = 0, S = 0.1, V = 0.1, initOpt = ThermoPower.Choices.Init.Options.steadyState, pstart = 11.2e5) annotation(
      Placement(transformation(extent = {{-38, -10}, {-18, 10}}, rotation = 0)));
    ThermoPower.Gas.SourceMassFlow Wfuel(redeclare package Medium = NaturalGas, use_in_w0 = true) annotation(
      Placement(transformation(extent = {{-50, 28}, {-30, 48}}, rotation = 0)));
    ThermoPower.Gas.PressDrop PressDrop1(redeclare package Medium = FlueGas, FFtype = ThermoPower.Choices.PressDrop.FFtypes.OpPoint, rhonom = 3.3, wnom = 158.9, pstart = 11.2e5, dpnom = 0.426e5) annotation(
      Placement(transformation(extent = {{-4, -10}, {16, 10}}, rotation = 0)));
    ThermoPower.Gas.SensT SensT1(redeclare package Medium = FlueGas) annotation(
      Placement(transformation(extent = {{26, -6}, {46, 14}}, rotation = 0)));
    Modelica.Blocks.Sources.Step Step1(startTime = 0.5, height = -0.3, offset = 3.1) annotation(
      Placement(transformation(extent = {{-78, 56}, {-58, 76}}, rotation = 0)));
    ThermoPower.Gas.ValveLin ValveLin1(redeclare package Medium = FlueGas, Kv = 161.1 / 9.77e5) annotation(
      Placement(transformation(extent = {{54, -10}, {74, 10}}, rotation = 0)));
    ThermoPower.Gas.SinkPressure SinkP1(redeclare package Medium = FlueGas) annotation(
      Placement(transformation(extent = {{84, -10}, {104, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Constant Constant1(k = 1) annotation(
      Placement(transformation(extent = {{22, 28}, {42, 48}}, rotation = 0)));
    inner ThermoPower.System system annotation(
      Placement(transformation(extent = {{80, 80}, {100, 100}})));
  equation
    connect(Wfuel.flange, CombustionChamber1.inf) annotation(
      Line(points = {{-30, 38}, {-28, 38}, {-28, 10}}, color = {159, 159, 223}, thickness = 0.5));
    connect(Wcompressor.flange, CombustionChamber1.ina) annotation(
      Line(points = {{-60, 0}, {-38, 0}}, color = {159, 159, 223}, thickness = 0.5));
    connect(CombustionChamber1.out, PressDrop1.inlet) annotation(
      Line(points = {{-18, 0}, {-4, 0}}, color = {159, 159, 223}, thickness = 0.5));
    connect(PressDrop1.outlet, SensT1.inlet) annotation(
      Line(points = {{16, 0}, {30, 0}}, color = {159, 159, 223}, thickness = 0.5));
    connect(Step1.y, Wfuel.in_w0) annotation(
      Line(points = {{-57, 66}, {-46, 66}, {-46, 43}}, color = {0, 0, 127}));
    connect(ValveLin1.outlet, SinkP1.flange) annotation(
      Line(points = {{74, 0}, {84, 0}}, color = {159, 159, 223}, thickness = 0.5));
    connect(SensT1.outlet, ValveLin1.inlet) annotation(
      Line(points = {{42, 0}, {54, 0}}, color = {159, 159, 223}, thickness = 0.5));
    connect(Constant1.y, ValveLin1.cmd) annotation(
      Line(points = {{43, 38}, {64, 38}, {64, 7}}, color = {0, 0, 127}));
    annotation(
      Documentation(info = "<html>
This model tests the <tt>CombustionChamber</tt> model. The model start at steady state. At time t = 0.5, the fuel flow rate is reduced by 10%.

<p>Simulate for 5s.
</html>"),
      experiment(StopTime = 5));
  end TestCC3;

  package FlueGas "flue gas"
    extends Modelica.Media.IdealGases.Common.MixtureGasNasa(mediumName = "FlueGas", data = {Modelica.Media.IdealGases.Common.SingleGasesData.O2, Modelica.Media.IdealGases.Common.SingleGasesData.Ar, Modelica.Media.IdealGases.Common.SingleGasesData.H2O, Modelica.Media.IdealGases.Common.SingleGasesData.CO2, Modelica.Media.IdealGases.Common.SingleGasesData.N2}, fluidConstants = {Modelica.Media.IdealGases.Common.FluidData.O2, Modelica.Media.IdealGases.Common.FluidData.Ar, Modelica.Media.IdealGases.Common.FluidData.H2O, Modelica.Media.IdealGases.Common.FluidData.CO2, Modelica.Media.IdealGases.Common.FluidData.N2}, substanceNames = {"Oxygen", "Argon", "Water", "Carbondioxide", "Nitrogen"}, reference_X = {0.23, 0.02, 0.01, 0.04, 0.7}, referenceChoice = Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C, excludeEnthalpyOfFormation = false);
  end FlueGas;

  model CombustionChamber2 "Combustion Chamber"
    extends CombustionChamberBase(redeclare package Air = MyTest2.Air "O2, H2O, Ar, N2", redeclare package Fuel = MyTest2.NaturalGas "N2, CO2, CH4", redeclare package Exhaust = MyTest2.FlueGas "O2, Ar, H2O, CO2, N2");
    Real wcomb(final quantity = "MolarFlowRate", unit = "mol/s") "Molar Combustion rate (CH4)";
    SI.PerUnit lambda "Stoichiometric ratio (>1 if air flow is greater than stoichiometric)";
  protected
    Air.MassFraction ina_X[Air.nXi] = inStream(ina.Xi_outflow);
    Fuel.MassFraction inf_X[Fuel.nXi] = inStream(inf.Xi_outflow);
  equation
    wcomb = inf.m_flow * inf_X[3] / Fuel.data[3].MM "Combustion molar flow rate";
    lambda = ina.m_flow * ina_X[1] / Air.data[1].MM / (2 * wcomb);
    assert(lambda >= 1, "Not enough oxygen flow");
    der(MX[1]) = ina.m_flow * ina_X[1] + out.m_flow * fluegas.X[1] - 2 * wcomb * Exhaust.data[1].MM "oxygen";
    der(MX[2]) = ina.m_flow * ina_X[3] + out.m_flow * fluegas.X[2] "argon";
    der(MX[3]) = ina.m_flow * ina_X[2] + out.m_flow * fluegas.X[3] + 2 * wcomb * Exhaust.data[3].MM "water";
    der(MX[4]) = inf.m_flow * inf_X[2] + out.m_flow * fluegas.X[4] + wcomb * Exhaust.data[4].MM "carbondioxide";
    der(MX[5]) = ina.m_flow * ina_X[4] + out.m_flow * fluegas.X[5] + inf.m_flow * inf_X[1] "nitrogen";
    annotation(
      Icon(graphics),
      Documentation(info = "<html>
This model extends the CombustionChamber Base model, with the definition of the gases.
<p>In particular, the air inlet uses the <tt>Media.Air</tt> medium model, the fuel input uses the <tt>Media.NaturalGas</tt> medium model, and the flue gas outlet uses the <tt>Medium.FlueGas</tt> medium model.
<p>The composition of the outgoing gas is determined by the mass balance of every component, taking into account the combustion reaction CH4+2O2--->2H2O+CO2.</p>
<p>The model assumes complete combustion, so that it is only valid if the oxygen flow at the air inlet is greater than the stoichiometric flow corresponding to the flow at the fuel inlet.</p>

</html>", revisions = "<html>
<ul>
<li><i>31 Jan 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
 Combustion Chamber model restructured using inheritance.
     <p>  First release.
 </li>
</ul>
</html>"));
  end CombustionChamber2;
end CombustionChamberTest2;