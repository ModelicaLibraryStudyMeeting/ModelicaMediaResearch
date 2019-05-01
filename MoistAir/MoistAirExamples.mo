model MoistAirExamples
  model SatPress
    replaceable package Medium = Modelica.Media.Air.MoistAir;
    Modelica.Blocks.Interfaces.RealInput T annotation(
      Placement(visible = true, transformation(origin = {-88, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-74, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    Modelica.Blocks.Interfaces.RealOutput saturationPressure1 annotation(
      Placement(visible = true, transformation(origin = {86, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {84, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Interfaces.RealOutput saturationPressure1_der annotation(
      Placement(visible = true, transformation(origin = {86, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {86, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    saturationPressure1 = Medium.saturationPressure(T + 273.15);
    saturationPressure1_der = Medium.saturationPressure_der(T + 273.15, 1.0);
    annotation(
      Diagram,
      Icon(graphics = {Rectangle(extent = {{-80, 80}, {80, -80}}), Text(origin = {-25, -16}, extent = {{-29, 28}, {-3, 6}}, textString = "T"), Text(origin = {136, 23}, extent = {{-100, 33}, {-74, 1}}, textString = "Psat"), Text(origin = {176, -119}, extent = {{-162, 91}, {-108, 67}}, textString = "Psat_der")}, coordinateSystem(initialScale = 0.1)));
  end SatPress;

  model SetPressTest
    MoistAirExamples.SatPress satPress1 annotation(
      Placement(visible = true, transformation(origin = {-1.11022e-15, 8.88178e-16}, extent = {{-26, -26}, {26, 26}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp temperature(duration = 455, height = 455, offset = -83, startTime = 0) annotation(
      Placement(visible = true, transformation(origin = {-60, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  equation
    connect(temperature.y, satPress1.T) annotation(
      Line(points = {{-38, 0}, {-19, 0}}, color = {0, 0, 127}));
    annotation(
      __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"),
      experiment(StartTime = 0, StopTime = 455, Tolerance = 1e-6, Interval = 1));
  end SetPressTest;

  model SatPressL
    replaceable package Medium = Modelica.Media.Air.MoistAir;
    Modelica.Blocks.Interfaces.RealInput T annotation(
      Placement(visible = true, transformation(origin = {-82, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-80, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    Modelica.Blocks.Interfaces.RealOutput p_sat annotation(
      Placement(visible = true, transformation(origin = {70, -8}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {82, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    p_sat = Medium.saturationPressureLiquid(T);
    annotation(
      Diagram,
      Icon(graphics = {Rectangle(extent = {{-80, 80}, {80, -80}}), Text(origin = {-35, -12}, extent = {{-29, 28}, {17, -12}}, textString = "T"), Text(origin = {42, -1}, extent = {{-36, 29}, {24, -23}}, textString = "p_sat")}));
  end SatPressL;

  model SetPressLTest
    SatPressL satPress1 annotation(
      Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp temperature(duration = 100, height = 100, offset = 273.16, startTime = 0) annotation(
      Placement(visible = true, transformation(origin = {-44, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(temperature.y, satPress1.T) annotation(
      Line(points = {{-32, 0}, {-8, 0}, {-8, 0}, {-8, 0}}, color = {0, 0, 127}));
    annotation(
      __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"));
  end SetPressLTest;

  model SatPressS
    replaceable package Medium = Modelica.Media.Air.MoistAir;
    Modelica.Blocks.Interfaces.RealInput T annotation(
      Placement(visible = true, transformation(origin = {-82, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-80, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    Modelica.Blocks.Interfaces.RealOutput p_sat annotation(
      Placement(visible = true, transformation(origin = {70, -8}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {82, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    p_sat = Medium.sublimationPressureIce(T);
    annotation(
      Diagram,
      Icon(graphics = {Rectangle(extent = {{-80, 80}, {80, -80}}), Text(origin = {-35, -12}, extent = {{-29, 28}, {17, -12}}, textString = "T"), Text(origin = {42, -1}, extent = {{-36, 29}, {24, -23}}, textString = "p_sat")}));
  end SatPressS;

  model SetPressSTest
    SatPressS satPress1 annotation(
      Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp temperature(duration = 73.15, height = 73.15, offset = 200, startTime = 0) annotation(
      Placement(visible = true, transformation(origin = {-44, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(temperature.y, satPress1.T) annotation(
      Line(points = {{-32, 0}, {-8, 0}, {-8, 0}, {-8, 0}}, color = {0, 0, 127}));
    annotation(
      __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"),
      experiment(StartTime = 0, StopTime = 73.15, Tolerance = 1e-6, Interval = 0.1463));
  end SetPressSTest;

  model VolumeCompression
    replaceable package Medium = MyMoistAir;
    Modelica.Fluid.Machines.SweptVolume sweptVolume1(redeclare package Medium = Medium, T_start = 293.15, X_start = {0.006, 0.994}, clearance = 0.1, energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyStateInitial, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, p_start = 101325, pistonCrossArea = 1, use_HeatTransfer = true, use_T_start = true) annotation(
      Placement(visible = true, transformation(origin = {-25, -33}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Mechanics.Translational.Sources.Position position1(v(fixed = true, start = 0)) annotation(
      Placement(visible = true, transformation(origin = {-24, 8}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Blocks.Sources.Ramp displacement(duration = 100, height = -0.9, offset = 0.9, startTime = 0) annotation(
      Placement(visible = true, transformation(origin = {-54, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature1(T = 293.15) annotation(
      Placement(visible = true, transformation(origin = {44, -32}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
    inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {30, 28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor1(G = 1e8) annotation(
      Placement(visible = true, transformation(origin = {10, -32}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
  equation
    connect(sweptVolume1.heatPort, thermalConductor1.port_b) annotation(
      Line(points = {{-15, -33}, {-7, -33}, {-7, -32}, {0, -32}}, color = {191, 0, 0}));
    connect(thermalConductor1.port_a, fixedTemperature1.port) annotation(
      Line(points = {{20, -32}, {34, -32}}, color = {191, 0, 0}));
    connect(position1.flange, sweptVolume1.flange) annotation(
      Line(points = {{-24, -2}, {-24, -20}, {-25, -20}, {-25, -23}}, color = {0, 127, 0}));
    connect(displacement.y, position1.s_ref) annotation(
      Line(points = {{-43, 38}, {-23, 38}, {-23, 17}, {-24, 17}, {-24, 20}}, color = {0, 0, 127}));
    annotation(
      __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"),
      experiment(StartTime = 0, StopTime = 100, Tolerance = 1e-6, Interval = 0.2));
  end VolumeCompression;

  package MyMoistAir
    extends Modelica.Media.Air.MoistAir;

    model extends BaseProperties
        MassFraction X_liquid1 "Mass fraction of liquid or solid water";
        MassFraction X_steam1 "Mass fraction of steam water";
        MassFraction X_air1 "Mass fraction of air";
        MassFraction X_sat1 "Steam water mass fraction of saturation boundary in kg_water/kg_moistair";
        MassFraction x_sat1 "Steam water mass content of saturation boundary in kg_water/kg_dryair";
        AbsolutePressure p_steam_sat1 "partial saturation pressure of steam";
        AbsolutePressure p_steam1;
        //            AbsolutePressure p_air1;
        //           Real p_diff;

      equation
        X_liquid1 = X_liquid;
        X_steam1 = X_steam;
        X_air1 = X_air;
        X_sat1 = X_sat;
        x_sat1 = x_sat;
        p_steam_sat1 = p_steam_sat;
      p_steam1 = d * X_steam * steam.R * T;
//            p_air1 = d * X_air * dryair.R * T;
//            p_diff = p - (p_air1 + p_steam1);
    end BaseProperties;
  equation

  end MyMoistAir;

      model HevT
              replaceable package Medium = Modelica.Media.Air.MoistAir;
              Modelica.Blocks.Interfaces.RealInput T annotation(
                      Placement(visible = true, transformation(origin = {-56, 2}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-74, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
              Modelica.Blocks.Interfaces.RealOutput hev annotation(
                      Placement(visible = true, transformation(origin = {48, 6}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {84, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
      hev = Medium.enthalpyOfVaporization(T+273.15);
              annotation(
                      Icon(graphics = {Rectangle(origin = {1, 1}, extent = {{-81, 79}, {79, -81}}), Text(origin = {-29, -16}, extent = {{-25, 36}, {15, -2}}, textString = "T"), Text(origin = {35, 5}, extent = {{-31, 67}, {31, -67}}, textString = "h_eva")}, coordinateSystem(initialScale = 0.1)),
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002),
  __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"));
      end HevT;
  
  

    model VolumeHeating
        replaceable package Medium = MyMoistAir;
        inner Modelica.Fluid.System system annotation(
            Placement(visible = true, transformation(origin = {74, 76}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Machines.SweptVolume sweptVolume1(redeclare package Medium = Medium, T_start = 268.15, X_start = {0.025, 0.975}, clearance = 0.1, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, nPorts = 0, p_start = 101325, pistonCrossArea = 1, use_HeatTransfer = true, use_T_start = true) annotation(
            Placement(visible = true, transformation(origin = {4, -24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Mechanics.Translational.Sensors.ForceSensor forceSensor1 annotation(
            Placement(visible = true, transformation(origin = {4, 24}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
        Modelica.Mechanics.Translational.Sources.Speed speed1(s(start = 0.9))  annotation(
            Placement(visible = true, transformation(origin = {4, 58}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
        Modelica.Blocks.Math.Gain gain1(k = 0.0002) annotation(
            Placement(visible = true, transformation(origin = {38, 54}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor1(G = 10)  annotation(
      Placement(visible = true, transformation(origin = {-28, -24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature1(T = 308.15)  annotation(
      Placement(visible = true, transformation(origin = {-68, -24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
    connect(gain1.y, speed1.v_ref) annotation(
      Line(points = {{38, 65}, {38, 88}, {4, 88}, {4, 70}}, color = {0, 0, 127}));
    connect(forceSensor1.f, gain1.u) annotation(
      Line(points = {{15, 16}, {37, 16}, {37, 42}}, color = {0, 0, 127}));
    connect(sweptVolume1.heatPort, thermalConductor1.port_b) annotation(
      Line(points = {{-6, -24}, {-18, -24}, {-18, -24}, {-18, -24}}, color = {191, 0, 0}));
    connect(forceSensor1.flange_a, sweptVolume1.flange) annotation(
      Line(points = {{4, 14}, {4, -14}}, color = {0, 127, 0}));
    connect(thermalConductor1.port_a, fixedTemperature1.port) annotation(
      Line(points = {{-38, -24}, {-58, -24}, {-58, -24}, {-58, -24}}, color = {191, 0, 0}));
    connect(speed1.flange, forceSensor1.flange_b) annotation(
      Line(points = {{4, 48}, {4, 34}}, color = {0, 127, 0}));
        annotation(
            __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"),
            experiment(StartTime = 0, StopTime = 2000, Tolerance = 1e-6, Interval = 0.1));
    end VolumeHeating;

  model HevTTest
  MoistAirExamples.HevT hevT1 annotation(
      Placement(visible = true, transformation(origin = {20, 0}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
  Modelica.Blocks.Sources.Ramp temperature(duration = 373, height = 373, offset = 0, startTime = 0)  annotation(
      Placement(visible = true, transformation(origin = {-16, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(temperature.y, hevT1.T) annotation(
      Line(points = {{-4, 0}, {8, 0}}, color = {0, 0, 127}));
  annotation(
      experiment(StartTime = 0, StopTime = 373, Tolerance = 1e-6, Interval = 0.746),
      __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"));end HevTTest;
equation

  annotation(
    uses(Modelica(version = "3.2.3")),
    __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"));
end MoistAirExamples;