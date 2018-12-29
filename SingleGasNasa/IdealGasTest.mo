package IdealGasTest
  import Modelica.Media.Interfaces.Choices.ReferenceEnthalpy;
  import SI = Modelica.SIunits;

  package EnthalpyTest
    model Check_h_T
      replaceable package Medium = Modelica.Media.IdealGases.SingleGases.CO2;
      parameter Boolean exclEnthForm = false;
      parameter ReferenceEnthalpy refChoice = ReferenceEnthalpy.ZeroAt25C;
      parameter Modelica.SIunits.SpecificEnthalpy h_off = 0.0;
      Modelica.Blocks.Interfaces.RealInput T annotation(
        Placement(visible = true, transformation(origin = {-68, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-68, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput h annotation(
        Placement(visible = true, transformation(origin = {80, 36}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {82, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput hmol annotation(
        Placement(visible = true, transformation(origin = {82, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {84, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      h = Modelica.Media.IdealGases.Common.Functions.h_T(Medium.data, T, exclEnthForm, refChoice, h_off);
      hmol = h * Medium.data.MM;
      annotation(
        Icon(graphics = {Rectangle(origin = {-3, 8}, extent = {{-77, 72}, {83, -88}}), Text(origin = {-42, 10}, extent = {{-2, 2}, {32, -32}}, textString = "T"), Text(origin = {19, -54}, extent = {{-11, -2}, {45, 30}}, textString = "hmol"), Text(origin = {37, 26}, extent = {{-11, -2}, {45, 30}}, textString = "h")}, coordinateSystem(initialScale = 0.1)));
    end Check_h_T;

    model Check_h_T_Test
      replaceable package Medium = Modelica.Media.IdealGases.SingleGases.CH4;
      Check_h_T check_h_T1(redeclare package Medium = Medium, exclEnthForm = false, refChoice = ReferenceEnthalpy.ZeroAt25C) annotation(
        Placement(visible = true, transformation(origin = {-18, 52}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Check_h_T check_h_T2(redeclare package Medium = Medium, exclEnthForm = true, refChoice = ReferenceEnthalpy.ZeroAt25C) annotation(
        Placement(visible = true, transformation(origin = {-18, 24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Check_h_T check_h_T3(redeclare package Medium = Medium, exclEnthForm = true, refChoice = ReferenceEnthalpy.ZeroAt0K) annotation(
        Placement(visible = true, transformation(origin = {-18, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Check_h_T check_h_T4(redeclare package Medium = Medium, exclEnthForm = true, h_off = 200000, refChoice = ReferenceEnthalpy.UserDefined) annotation(
        Placement(visible = true, transformation(origin = {-18, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Ramp temperature(duration = 1000, height = 800, offset = 200, startTime = 0) annotation(
        Placement(visible = true, transformation(origin = {-60, 8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(temperature.y, check_h_T4.T) annotation(
        Line(points = {{-48, 8}, {-40, 8}, {-40, -30}, {-24, -30}, {-24, -30}}, color = {0, 0, 127}));
      connect(temperature.y, check_h_T3.T) annotation(
        Line(points = {{-48, 8}, {-40, 8}, {-40, -2}, {-24, -2}, {-24, -2}}, color = {0, 0, 127}));
      connect(temperature.y, check_h_T2.T) annotation(
        Line(points = {{-48, 8}, {-40, 8}, {-40, 24}, {-24, 24}, {-24, 24}}, color = {0, 0, 127}));
      connect(temperature.y, check_h_T1.T) annotation(
        Line(points = {{-48, 8}, {-40, 8}, {-40, 52}, {-24, 52}, {-24, 52}}, color = {0, 0, 127}));
      annotation(
        experiment(StartTime = 0, StopTime = 1000, Tolerance = 1e-6, Interval = 2),
        __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"));
    end Check_h_T_Test;

    model Check_h_T_Test2
      replaceable package Medium = Modelica.Media.IdealGases.SingleGases.CH4;
      Check_h_T check_h_T1(redeclare package Medium = Medium, exclEnthForm = false, refChoice = ReferenceEnthalpy.ZeroAt25C) annotation(
        Placement(visible = true, transformation(origin = {-18, 52}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Check_h_T check_h_T2(redeclare package Medium = Medium, exclEnthForm = true, refChoice = ReferenceEnthalpy.ZeroAt25C) annotation(
        Placement(visible = true, transformation(origin = {-18, 24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Check_h_T check_h_T3(redeclare package Medium = Medium, exclEnthForm = true, refChoice = ReferenceEnthalpy.ZeroAt0K) annotation(
        Placement(visible = true, transformation(origin = {-18, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Check_h_T check_h_T4(redeclare package Medium = Medium, exclEnthForm = true, h_off = 200000, refChoice = ReferenceEnthalpy.UserDefined) annotation(
        Placement(visible = true, transformation(origin = {-18, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant const(k = 300)  annotation(
        Placement(visible = true, transformation(origin = {-62, 22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(const.y, check_h_T4.T) annotation(
        Line(points = {{-50, 22}, {-40, 22}, {-40, -30}, {-24, -30}, {-24, -30}}, color = {0, 0, 127}));
      connect(const.y, check_h_T3.T) annotation(
        Line(points = {{-50, 22}, {-40, 22}, {-40, -2}, {-24, -2}, {-24, -2}}, color = {0, 0, 127}));
      connect(const.y, check_h_T2.T) annotation(
        Line(points = {{-50, 22}, {-40, 22}, {-40, 24}, {-24, 24}, {-24, 24}}, color = {0, 0, 127}));
      connect(const.y, check_h_T1.T) annotation(
        Line(points = {{-50, 22}, {-40, 22}, {-40, 52}, {-24, 52}, {-24, 52}}, color = {0, 0, 127}));
      annotation(
        experiment(StartTime = 0, StopTime = 1000, Tolerance = 1e-6, Interval = 2),
        __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"));
    end Check_h_T_Test2;
  end EnthalpyTest;

  package StateSelectTest
        model Room
            replaceable package Medium = Modelica.Media.Air.DryAirNasa;
            parameter SI.Volume V = 22.0;
            Modelica.Fluid.Interfaces.FluidPort_b port_a(redeclare package Medium = Medium) annotation(
                Placement(visible = true, transformation(origin = {-46, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-60, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
            Modelica.Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium = Medium) annotation(
                Placement(visible = true, transformation(origin = {48, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {60, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
            Medium.BaseProperties medium;
            SI.Mass M;
            SI.Energy U;
        equation
            M = medium.d * V;
            U = medium.u * M;
            der(M) = port_a.m_flow + port_b.m_flow;
            der(U) = actualStream(port_a.h_outflow) * port_a.m_flow + actualStream(port_b.h_outflow) * port_b.m_flow;
            port_a.p = medium.p;
            port_b.p = medium.p;
            port_a.h_outflow = medium.h;
            port_b.h_outflow = medium.h;
        initial equation
            medium.T = 293.15;
            medium.p = 101325;
            annotation(
                Icon(graphics = {Rectangle(origin = {-20, 4}, extent = {{-40, 56}, {80, -64}}), Text(origin = {4, 75}, extent = {{-70, 15}, {54, -11}}, textString = "%name"), Polygon(origin = {-1.02, 49}, fillColor = {255, 128, 0}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, points = {{-58.983, 11}, {-78.983, -21}, {81.017, -21}, {61.017, 11}, {61.017, 11}, {-0.982991, 11}, {-58.983, 11}})}, coordinateSystem(initialScale = 0.1)),
                experiment(StartTime = 0, StopTime = 5000, Tolerance = 1e-06, Interval = 10),
                __OpenModelica_simulationFlags(lv = "LOG_LS,LOG_NLS,LOG_STATS", outputFormat = "mat", s = "dassl"));
        end Room;
    

    model RoomTest_prefer
      replaceable package Medium = Modelica.Media.Air.DryAirNasa(AbsolutePressure(nominal = 100000.0));
      IdealGasTest.StateSelectTest.Room room1(redeclare package Medium = Medium, medium(preferredMediumStates = true)) annotation(
        Placement(visible = true, transformation(origin = {0, 2}, extent = {{-18, -18}, {18, 18}}, rotation = 0)));
      Modelica.Fluid.Sources.MassFlowSource_T boundary1(redeclare package Medium = Medium, T = 323.15, m_flow = 0.05, nPorts = 1) annotation(
        Placement(visible = true, transformation(origin = {-46, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Fluid.Sources.MassFlowSource_T boundary2(redeclare package Medium = Medium, m_flow = -0.045, nPorts = 1) annotation(
        Placement(visible = true, transformation(origin = {50, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    equation
      connect(boundary2.ports[1], room1.port_b) annotation(
        Line(points = {{40, 0}, {12, 0}, {12, 0}, {10, 0}}, color = {0, 127, 255}));
      connect(boundary1.ports[1], room1.port_a) annotation(
        Line(points = {{-36, 0}, {-11, 0}}, color = {0, 127, 255}));
      annotation(
        experiment(StartTime = 0, StopTime = 5000, Tolerance = 1e-06, Interval = 10),
        __OpenModelica_simulationFlags(outputFormat = "mat", s = "dassl", lv = "LOG_LS,LOG_NLS,LOG_STATS"));
    end RoomTest_prefer;


    model RoomTest_default
      replaceable package Medium = Modelica.Media.Air.DryAirNasa(AbsolutePressure(nominal = 100000.0));
      IdealGasTest.StateSelectTest.Room room1(redeclare package Medium = Medium, medium(preferredMediumStates = false)) annotation(
        Placement(visible = true, transformation(origin = {0, 2}, extent = {{-18, -18}, {18, 18}}, rotation = 0)));
      Modelica.Fluid.Sources.MassFlowSource_T boundary1(redeclare package Medium = Medium, T = 323.15, m_flow = 0.05, nPorts = 1) annotation(
        Placement(visible = true, transformation(origin = {-46, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Fluid.Sources.MassFlowSource_T boundary2(redeclare package Medium = Medium, m_flow = -0.045, nPorts = 1) annotation(
        Placement(visible = true, transformation(origin = {50, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    equation
      connect(boundary2.ports[1], room1.port_b) annotation(
        Line(points = {{40, 0}, {12, 0}, {12, 0}, {10, 0}}, color = {0, 127, 255}));
      connect(boundary1.ports[1], room1.port_a) annotation(
        Line(points = {{-36, 0}, {-11, 0}}, color = {0, 127, 255}));
      annotation(
        experiment(StartTime = 0, StopTime = 5000, Tolerance = 1e-06, Interval = 10),
        __OpenModelica_simulationFlags(outputFormat = "mat", s = "dassl", lv = "LOG_LS,LOG_NLS,LOG_STATS"));
    end RoomTest_default;

  end StateSelectTest;

  package EntropyTest
    model Check_s0_T
      replaceable package Medium = Modelica.Media.IdealGases.SingleGases.CO2;
      Medium.ThermodynamicState state;
      Modelica.Blocks.Interfaces.RealInput T annotation(
        Placement(visible = true, transformation(origin = {-64, -30}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-66, -36}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput s0 annotation(
        Placement(visible = true, transformation(origin = {72, 56}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {84, 46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Interfaces.RealInput p annotation(
        Placement(visible = true, transformation(origin = {-66, 48}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-66, 40}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    Modelica.Blocks.Interfaces.RealOutput s annotation(
        Placement(visible = true, transformation(origin = {74, 6}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {84, -8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput smol annotation(
        Placement(visible = true, transformation(origin = {76, -44}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {82, -54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      state = Medium.setState_pT(p, T);
      s = Medium.specificEntropy(state);
      smol = s*Medium.molarMass(state);
      s0 = Modelica.Media.IdealGases.Common.Functions.s0_T(Medium.data, T);
      annotation(
        Icon(graphics = {Rectangle(origin = {-3, 8}, extent = {{-77, 72}, {83, -88}}), Text(origin = {-40, -24}, extent = {{-2, 2}, {32, -32}}, textString = "T"),  Text(origin = {25, 28}, extent = {{-11, -2}, {45, 30}}, textString = "s0"), Text(origin = {31, -22}, extent = {{-11, -2}, {45, 30}}, textString = "s"), Text(origin = {-16, 26}, extent = {{-24, 32}, {10, -4}}, textString = "p"), Text(origin = {15, -62}, extent = {{-11, -2}, {45, 30}}, textString = "smol")}, coordinateSystem(initialScale = 0.1)));
    end Check_s0_T;



    model Check_s0_T_Test
      replaceable package Medium = Modelica.Media.IdealGases.SingleGases.H2O;
      Check_s0_T check_s0_T(redeclare package Medium = Medium) annotation(
        Placement(visible = true, transformation(origin = {-16, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Ramp temperature(duration = 1000, height = 800.0, offset = 200.0, startTime = 0) annotation(
        Placement(visible = true, transformation(origin = {-62, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant pressure1(k = 101325)  annotation(
        Placement(visible = true, transformation(origin = {-64, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant pressure2(k = 1001325) annotation(
        Placement(visible = true, transformation(origin = {12, 44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Check_s0_T check_s0_T1(redeclare package Medium = Medium) annotation(
        Placement(visible = true, transformation(origin = {60, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(temperature.y, check_s0_T1.T) annotation(
        Line(points = {{-50, 2}, {34, 2}, {34, 14}, {54, 14}, {54, 14}}, color = {0, 0, 127}));
      connect(check_s0_T1.p, pressure2.y) annotation(
        Line(points = {{53, 22}, {34, 22}, {34, 44}, {24, 44}}, color = {0, 0, 127}));
      connect(temperature.y, check_s0_T.T) annotation(
        Line(points = {{-51, 2}, {-32, 2}, {-32, 14}, {-23, 14}}, color = {0, 0, 127}));
      connect(pressure1.y, check_s0_T.p) annotation(
        Line(points = {{-53, 42}, {-32, 42}, {-32, 22}, {-23, 22}}, color = {0, 0, 127}));
      annotation(
        experiment(StartTime = 0, StopTime = 1000, Tolerance = 1e-6, Interval = 2),
        __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"));
    end Check_s0_T_Test;














  end EntropyTest;

  package TemperatureEstimation
    model T_ph_Check
      replaceable package Medium = Modelica.Media.IdealGases.SingleGases.CH4;
      Medium.ThermodynamicState state;
      Modelica.Blocks.Interfaces.RealInput p annotation(
        Placement(visible = true, transformation(origin = {-66, 40}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-76, 40}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput h annotation(
        Placement(visible = true, transformation(origin = {-66, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-76, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput T annotation(
        Placement(visible = true, transformation(origin = {66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {80, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      state = Medium.setState_ph(p, h);
      T = state.T;
      annotation(
        Diagram,
        Icon(graphics = {Rectangle(origin = {0, 1}, extent = {{-80, 79}, {80, -81}}), Text(origin = {-31, 43}, extent = {{-19, 23}, {19, -23}}, textString = "p"), Text(origin = {-31, -43}, extent = {{-19, 23}, {19, -23}}, textString = "h"), Text(origin = {41, -1}, extent = {{-19, 23}, {19, -23}}, textString = "T")}));
    end T_ph_Check;

    model H_pT_Check
      replaceable package Medium = Modelica.Media.IdealGases.SingleGases.CH4;
      Medium.ThermodynamicState state;
      Modelica.Blocks.Interfaces.RealInput p annotation(
        Placement(visible = true, transformation(origin = {-66, 40}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-76, 40}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput T annotation(
        Placement(visible = true, transformation(origin = {-66, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-76, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput h annotation(
        Placement(visible = true, transformation(origin = {66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {80, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      state = Medium.setState_pT(p, T);
      h = Medium.specificEnthalpy(state);
      annotation(
        Diagram,
        Icon(graphics = {Rectangle(origin = {0, 1}, extent = {{-80, 79}, {80, -81}}), Text(origin = {-31, 43}, extent = {{-19, 23}, {19, -23}}, textString = "p"), Text(origin = {51, -1}, extent = {{-19, 23}, {19, -23}}, textString = "h"), Text(origin = {-35, -41}, extent = {{-19, 23}, {19, -23}}, textString = "T")}, coordinateSystem(initialScale = 0.1)));
    end H_pT_Check;


    model T_ph_Test
      H_pT_Check h_pT_Check1 annotation(
        Placement(visible = true, transformation(origin = {2, 44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant pressure(k = 101325) annotation(
        Placement(visible = true, transformation(origin = {-46, 54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Ramp temperature(duration = 1000, height = 800, offset = 200) annotation(
        Placement(visible = true, transformation(origin = {-46, 22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      T_ph_Check t_ph_Check1 annotation(
        Placement(visible = true, transformation(origin = {62, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(h_pT_Check1.h, t_ph_Check1.h) annotation(
        Line(points = {{10, 44}, {36, 44}, {36, 46}, {54, 46}}, color = {0, 0, 127}));
      connect(pressure.y, t_ph_Check1.p) annotation(
        Line(points = {{-34, 54}, {54, 54}}, color = {0, 0, 127}));
      connect(temperature.y, h_pT_Check1.T) annotation(
        Line(points = {{-34, 22}, {-20, 22}, {-20, 40}, {-6, 40}, {-6, 40}}, color = {0, 0, 127}));
      connect(pressure.y, h_pT_Check1.p) annotation(
        Line(points = {{-34, 54}, {-20, 54}, {-20, 48}, {-6, 48}, {-6, 48}}, color = {0, 0, 127}));
      annotation(
        __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"));
    end T_ph_Test;



    model S_pT_Check
      replaceable package Medium = Modelica.Media.IdealGases.SingleGases.CH4;
      Medium.ThermodynamicState state;
      Modelica.Blocks.Interfaces.RealInput p annotation(
        Placement(visible = true, transformation(origin = {-66, 40}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-76, 40}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput T annotation(
        Placement(visible = true, transformation(origin = {-66, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-76, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput s annotation(
        Placement(visible = true, transformation(origin = {66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {80, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      state = Medium.setState_pT(p, T);
      s = Medium.specificEntropy(state);
      annotation(
        Diagram,
        Icon(graphics = {Rectangle(origin = {0, 1}, extent = {{-80, 79}, {80, -81}}), Text(origin = {-31, 43}, extent = {{-19, 23}, {19, -23}}, textString = "p"), Text(origin = {51, -1}, extent = {{-19, 23}, {19, -23}}, textString = "s"), Text(origin = {-35, -41}, extent = {{-19, 23}, {19, -23}}, textString = "T")}, coordinateSystem(initialScale = 0.1)));
    end S_pT_Check;


    model T_ps_Check
      replaceable package Medium = Modelica.Media.IdealGases.SingleGases.CH4;
      Medium.ThermodynamicState state;
      Modelica.Blocks.Interfaces.RealInput p annotation(
        Placement(visible = true, transformation(origin = {-66, 40}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-76, 40}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput s annotation(
        Placement(visible = true, transformation(origin = {-66, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-76, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput T annotation(
        Placement(visible = true, transformation(origin = {66, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {80, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      state = Medium.setState_ps(p, s);
      T = state.T;
      annotation(
        Diagram,
        Icon(graphics = {Rectangle(origin = {0, 1}, extent = {{-80, 79}, {80, -81}}), Text(origin = {-31, 43}, extent = {{-19, 23}, {19, -23}}, textString = "p"), Text(origin = {-31, -43}, extent = {{-19, 23}, {19, -23}}, textString = "s"), Text(origin = {41, -1}, extent = {{-19, 23}, {19, -23}}, textString = "T")}, coordinateSystem(initialScale = 0.1)));
    end T_ps_Check;


    model T_ps_Test
      S_pT_Check s_pT_Check1 annotation(
        Placement(visible = true, transformation(origin = {2, 44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant pressure(k = 101325) annotation(
        Placement(visible = true, transformation(origin = {-46, 54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Ramp temperature(duration = 1000, height = 800, offset = 200) annotation(
        Placement(visible = true, transformation(origin = {-46, 22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      T_ps_Check t_ps_Check1 annotation(
        Placement(visible = true, transformation(origin = {62, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(s_pT_Check1.s, t_ps_Check1.s) annotation(
        Line(points = {{10, 44}, {36, 44}, {36, 46}, {54, 46}}, color = {0, 0, 127}));
      connect(pressure.y, t_ps_Check1.p) annotation(
        Line(points = {{-34, 54}, {54, 54}}, color = {0, 0, 127}));
      connect(temperature.y, s_pT_Check1.T) annotation(
        Line(points = {{-34, 22}, {-20, 22}, {-20, 40}, {-6, 40}, {-6, 40}}, color = {0, 0, 127}));
      connect(pressure.y, s_pT_Check1.p) annotation(
        Line(points = {{-34, 54}, {-20, 54}, {-20, 48}, {-6, 48}, {-6, 48}}, color = {0, 0, 127}));
      annotation(
        __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"));
    end T_ps_Test;








  end TemperatureEstimation;

  package EtaLamba
    model ELC_Check
      replaceable package Medium = Modelica.Media.IdealGases.SingleGases.CH4;
      Medium.ThermodynamicState state;
      Modelica.Blocks.Interfaces.RealInput p annotation(
        Placement(visible = true, transformation(origin = {-68, 40}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-76, 40}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput T annotation(
        Placement(visible = true, transformation(origin = {-66, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-74, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput eta annotation(
        Placement(visible = true, transformation(origin = {56, 54}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {82, 44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput lambda annotation(
        Placement(visible = true, transformation(origin = {40, 8}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {84, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput cp annotation(
        Placement(visible = true, transformation(origin = {54, -74}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {84, -58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput lambda2 annotation(
        Placement(visible = true, transformation(origin = {82, -24}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {84, -8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput lambda3 annotation(
        Placement(visible = true, transformation(origin = {80, -52}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {84, -34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput y annotation(
        Placement(visible = true, transformation(origin = {86, 86}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {78, 68}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      state = Medium.setState_pT(p, T);
      eta = Medium.dynamicViscosity(state);
      lambda = Medium.thermalConductivity(state = state, method = 2);
      lambda2 = thermalConductivityEstimate2(Medium.specificHeatCapacityCp(state), eta, method = 2, data = Medium.data);
      lambda3 = thermalConductivityEstimate3(Medium.specificHeatCapacityCp(state), eta, method = 2, data = Medium.data);
      cp = Medium.specificHeatCapacityCp(state) * Medium.molarMass(state);
      y = Medium.heatCapacity_cp(state) * Medium.data.MM / Modelica.Constants.R - Medium.heatCapacity_cp(state) / Medium.data.R;
      annotation(
        Icon(graphics = {Rectangle(extent = {{-80, 80}, {80, -80}}), Text(origin = {-35, 38}, extent = {{-23, 22}, {23, -22}}, textString = "p"), Text(origin = {-35, -36}, extent = {{-23, 22}, {23, -22}}, textString = "T"), Text(origin = {37, 54}, extent = {{-11, 2}, {23, -22}}, textString = "eta"), Text(origin = {43, 0}, extent = {{-55, 50}, {23, -22}}, textString = "lambda"), Text(origin = {47, -50}, extent = {{-11, 2}, {23, -22}}, textString = "cp"), Text(origin = {43, -24}, extent = {{-55, 50}, {23, -22}}, textString = "lambda2"), Text(origin = {43, -50}, extent = {{-55, 50}, {23, -22}}, textString = "lambda3")}, coordinateSystem(initialScale = 0.1)),
        __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"));end ELC_Check;














    model ETC_Test
    IdealGasTest.EtaLamba.ELC_Check eLC_Check1 annotation(
        Placement(visible = true, transformation(origin = {4, 62}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant pressure(k = 101325)  annotation(
        Placement(visible = true, transformation(origin = {-50, 86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant temperature(k = 150) annotation(
        Placement(visible = true, transformation(origin = {-48, 58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  IdealGasTest.EtaLamba.ELC_Check eLC_Check2 annotation(
        Placement(visible = true, transformation(origin = {4, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant const(k = 200) annotation(
        Placement(visible = true, transformation(origin = {-48, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant const1(k = 250) annotation(
        Placement(visible = true, transformation(origin = {-48, -6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  IdealGasTest.EtaLamba.ELC_Check eLC_Check3 annotation(
        Placement(visible = true, transformation(origin = {4, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant const2(k = 300) annotation(
        Placement(visible = true, transformation(origin = {-48, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  IdealGasTest.EtaLamba.ELC_Check eLC_Check4 annotation(
        Placement(visible = true, transformation(origin = {4, -36}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant const3(k = 350) annotation(
        Placement(visible = true, transformation(origin = {-48, -76}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  IdealGasTest.EtaLamba.ELC_Check eLC_Check5 annotation(
        Placement(visible = true, transformation(origin = {4, -72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant const4(k = 400) annotation(
        Placement(visible = true, transformation(origin = {52, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant const5(k = 450) annotation(
        Placement(visible = true, transformation(origin = {52, 28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  IdealGasTest.EtaLamba.ELC_Check eLC_Check6 annotation(
        Placement(visible = true, transformation(origin = {104, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  IdealGasTest.EtaLamba.ELC_Check eLC_Check7 annotation(
        Placement(visible = true, transformation(origin = {104, 32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  IdealGasTest.EtaLamba.ELC_Check eLC_Check8 annotation(
        Placement(visible = true, transformation(origin = {104, 64}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant const6(k = 500) annotation(
        Placement(visible = true, transformation(origin = {52, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  IdealGasTest.EtaLamba.ELC_Check eLC_Check9 annotation(
        Placement(visible = true, transformation(origin = {104, -34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant const7(k = 550) annotation(
        Placement(visible = true, transformation(origin = {52,-38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  IdealGasTest.EtaLamba.ELC_Check eLC_Check10 annotation(
        Placement(visible = true, transformation(origin = {104, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant const8(k = 600) annotation(
        Placement(visible = true, transformation(origin = {52, -74}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(pressure.y, eLC_Check10.p) annotation(
        Line(points = {{-38, 86}, {80, 86}, {80, -66}, {96, -66}, {96, -66}}, color = {0, 0, 127}));
      connect(pressure.y, eLC_Check9.p) annotation(
        Line(points = {{-38, 86}, {80, 86}, {80, -30}, {96, -30}, {96, -30}}, color = {0, 0, 127}));
      connect(pressure.y, eLC_Check6.p) annotation(
        Line(points = {{-38, 86}, {80, 86}, {80, 4}, {96, 4}, {96, 4}}, color = {0, 0, 127}));
      connect(pressure.y, eLC_Check7.p) annotation(
        Line(points = {{-38, 86}, {80, 86}, {80, 36}, {96, 36}, {96, 36}}, color = {0, 0, 127}));
      connect(pressure.y, eLC_Check8.p) annotation(
        Line(points = {{-38, 86}, {80, 86}, {80, 68}, {96, 68}, {96, 68}}, color = {0, 0, 127}));
      connect(pressure.y, eLC_Check5.p) annotation(
        Line(points = {{-38, 86}, {-20, 86}, {-20, -68}, {-4, -68}, {-4, -68}}, color = {0, 0, 127}));
      connect(pressure.y, eLC_Check4.p) annotation(
        Line(points = {{-38, 86}, {-20, 86}, {-20, -32}, {-4, -32}, {-4, -32}}, color = {0, 0, 127}));
      connect(pressure.y, eLC_Check3.p) annotation(
        Line(points = {{-38, 86}, {-20, 86}, {-20, 2}, {-4, 2}, {-4, 2}}, color = {0, 0, 127}));
      connect(pressure.y, eLC_Check2.p) annotation(
        Line(points = {{-38, 86}, {-20, 86}, {-20, 34}, {-4, 34}, {-4, 34}}, color = {0, 0, 127}));
      connect(pressure.y, eLC_Check1.p) annotation(
        Line(points = {{-38, 86}, {-20, 86}, {-20, 66}, {-4, 66}, {-4, 66}, {-4, 66}}, color = {0, 0, 127}));
      connect(const8.y, eLC_Check10.T) annotation(
        Line(points = {{64, -74}, {96, -74}, {96, -74}, {96, -74}}, color = {0, 0, 127}));
      connect(const7.y, eLC_Check9.T) annotation(
        Line(points = {{64, -38}, {96, -38}, {96, -38}, {96, -38}}, color = {0, 0, 127}));
      connect(const6.y, eLC_Check6.T) annotation(
        Line(points = {{64, -4}, {96, -4}, {96, -4}, {96, -4}}, color = {0, 0, 127}));
      connect(const5.y, eLC_Check7.T) annotation(
        Line(points = {{64, 28}, {96, 28}, {96, 28}, {96, 28}}, color = {0, 0, 127}));
      connect(const4.y, eLC_Check8.T) annotation(
        Line(points = {{64, 60}, {96, 60}, {96, 60}, {96, 60}}, color = {0, 0, 127}));
      connect(const3.y, eLC_Check5.T) annotation(
        Line(points = {{-36, -76}, {-4, -76}, {-4, -76}, {-4, -76}}, color = {0, 0, 127}));
      connect(temperature.y, eLC_Check1.T) annotation(
        Line(points = {{-37, 58}, {-3, 58}}, color = {0, 0, 127}));
      connect(const.y, eLC_Check2.T) annotation(
        Line(points = {{-37, 26}, {-3, 26}}, color = {0, 0, 127}));
      connect(const1.y, eLC_Check3.T) annotation(
        Line(points = {{-37, -6}, {-3, -6}}, color = {0, 0, 127}));
    connect(const2.y, eLC_Check4.T) annotation(
        Line(points = {{-37, -40}, {-3, -40}}, color = {0, 0, 127}));
    annotation(
        __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"));end ETC_Test;

    function thermalConductivityEstimate2 "Thermal conductivity of polyatomic gases(Eucken and Modified Eucken correlation)"
      extends Modelica.Icons.Function;
      input Modelica.Media.Interfaces.Types.SpecificHeatCapacity Cp "Constant pressure heat capacity";
      input Modelica.Media.Interfaces.Types.DynamicViscosity eta "Dynamic viscosity";
      input Integer method(min = 1, max = 2) = 1 "1: Eucken Method, 2: Modified Eucken Method";
      input Modelica.Media.IdealGases.Common.DataRecord data "Ideal gas data";
      output Modelica.Media.Interfaces.Types.ThermalConductivity lambda "Thermal conductivity [W/(m.k)]";
    algorithm
      lambda := if method == 1 then eta * (Cp - data.R + 9 / 4 * data.R) else eta * (Cp - data.R) * (1.32 + 1.77 / (Cp*data.MM/ Modelica.Constants.R- 1.0));
      annotation(
        smoothOrder = 2,
        Documentation(info = "<html>
<p>
This function provides two similar methods for estimating the
thermal conductivity of polyatomic gases.
The Eucken method (input method == 1) gives good results for low temperatures,
but it tends to give an underestimated value of the thermal conductivity
(lambda) at higher temperatures.<br>
The Modified Eucken method (input method == 2) gives good results for
high-temperatures, but it tends to give an overestimated value of the
thermal conductivity (lambda) at low temperatures.
</p>
</html>"));
    end thermalConductivityEstimate2;

    function thermalConductivityEstimate3 "Thermal conductivity of polyatomic gases(Eucken and Modified Eucken correlation)"
      extends Modelica.Icons.Function;
      input Modelica.Media.Interfaces.Types.SpecificHeatCapacity Cp "Constant pressure heat capacity";
      input Modelica.Media.Interfaces.Types.DynamicViscosity eta "Dynamic viscosity";
      input Integer method(min = 1, max = 2) = 1 "1: Eucken Method, 2: Modified Eucken Method";
      input Modelica.Media.IdealGases.Common.DataRecord data "Ideal gas data";
      output Modelica.Media.Interfaces.Types.ThermalConductivity lambda "Thermal conductivity [W/(m.k)]";
    algorithm
      lambda := if method == 1 then eta * (Cp - data.R + 9 / 4 * data.R) else eta * (Cp - data.R) * (1.32 + 1.77 / (Cp / data.R - 1.0));
      annotation(
        smoothOrder = 2,
        Documentation(info = "<html>
<p>
This function provides two similar methods for estimating the
thermal conductivity of polyatomic gases.
The Eucken method (input method == 1) gives good results for low temperatures,
but it tends to give an underestimated value of the thermal conductivity
(lambda) at higher temperatures.<br>
The Modified Eucken method (input method == 2) gives good results for
high-temperatures, but it tends to give an overestimated value of the
thermal conductivity (lambda) at low temperatures.
</p>
</html>"));
    end thermalConductivityEstimate3;

    model ETC_Test2
      IdealGasTest.EtaLamba.ELC_Check eLC_Check1(redeclare package Medium = Modelica.Media.IdealGases.SingleGases.Ar)  annotation(
        Placement(visible = true, transformation(origin = {4, 62}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant pressure(k = 101325) annotation(
        Placement(visible = true, transformation(origin = {-50, 86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      IdealGasTest.EtaLamba.ELC_Check eLC_Check2(redeclare package Medium = Modelica.Media.IdealGases.SingleGases.CH4)  annotation(
        Placement(visible = true, transformation(origin = {4, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      IdealGasTest.EtaLamba.ELC_Check eLC_Check3(redeclare package Medium = Modelica.Media.IdealGases.SingleGases.CH3OH)  annotation(
        Placement(visible = true, transformation(origin = {4, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      IdealGasTest.EtaLamba.ELC_Check eLC_Check4(redeclare package Medium = Modelica.Media.IdealGases.SingleGases.CO)  annotation(
        Placement(visible = true, transformation(origin = {4, -36}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      IdealGasTest.EtaLamba.ELC_Check eLC_Check5(redeclare package Medium = Modelica.Media.IdealGases.SingleGases.CO2)  annotation(
        Placement(visible = true, transformation(origin = {4, -72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      IdealGasTest.EtaLamba.ELC_Check eLC_Check6(redeclare package Medium = Modelica.Media.IdealGases.SingleGases.C2H2_vinylidene)  annotation(
        Placement(visible = true, transformation(origin = {104, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      IdealGasTest.EtaLamba.ELC_Check eLC_Check7(redeclare package Medium = Modelica.Media.IdealGases.SingleGases.C2H5OH)  annotation(
        Placement(visible = true, transformation(origin = {104, 32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      IdealGasTest.EtaLamba.ELC_Check eLC_Check8(redeclare package Medium = Modelica.Media.IdealGases.SingleGases.C2H6)  annotation(
        Placement(visible = true, transformation(origin = {104, 64}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      IdealGasTest.EtaLamba.ELC_Check eLC_Check9(redeclare package Medium = Modelica.Media.IdealGases.SingleGases.C3H6_propylene)  annotation(
        Placement(visible = true, transformation(origin = {104, -34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      IdealGasTest.EtaLamba.ELC_Check eLC_Check10(redeclare package Medium = Modelica.Media.IdealGases.SingleGases.C3H8)  annotation(
        Placement(visible = true, transformation(origin = {104, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Ramp ramp1(duration = 5800, height = 5800, offset = 200, startTime = 0)  annotation(
        Placement(visible = true, transformation(origin = {-50, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(ramp1.y, eLC_Check10.T) annotation(
        Line(points = {{-38, -80}, {66, -80}, {66, -74}, {96, -74}, {96, -74}}, color = {0, 0, 127}));
      connect(ramp1.y, eLC_Check9.T) annotation(
        Line(points = {{-38, -80}, {66, -80}, {66, -38}, {96, -38}, {96, -38}}, color = {0, 0, 127}));
      connect(ramp1.y, eLC_Check6.T) annotation(
        Line(points = {{-38, -80}, {66, -80}, {66, -4}, {96, -4}, {96, -4}}, color = {0, 0, 127}));
      connect(ramp1.y, eLC_Check7.T) annotation(
        Line(points = {{-38, -80}, {66, -80}, {66, 28}, {96, 28}, {96, 28}}, color = {0, 0, 127}));
      connect(ramp1.y, eLC_Check8.T) annotation(
        Line(points = {{-38, -80}, {66, -80}, {66, 60}, {96, 60}, {96, 60}}, color = {0, 0, 127}));
      connect(ramp1.y, eLC_Check5.T) annotation(
        Line(points = {{-38, -80}, {-30, -80}, {-30, -76}, {-4, -76}, {-4, -76}}, color = {0, 0, 127}));
      connect(ramp1.y, eLC_Check4.T) annotation(
        Line(points = {{-38, -80}, {-30, -80}, {-30, -40}, {-4, -40}, {-4, -40}}, color = {0, 0, 127}));
      connect(ramp1.y, eLC_Check3.T) annotation(
        Line(points = {{-38, -80}, {-30, -80}, {-30, -6}, {-4, -6}, {-4, -6}}, color = {0, 0, 127}));
      connect(ramp1.y, eLC_Check2.T) annotation(
        Line(points = {{-38, -80}, {-30, -80}, {-30, 26}, {-4, 26}, {-4, 26}}, color = {0, 0, 127}));
      connect(ramp1.y, eLC_Check1.T) annotation(
        Line(points = {{-38, -80}, {-30, -80}, {-30, 58}, {-4, 58}, {-4, 58}}, color = {0, 0, 127}));
      connect(pressure.y, eLC_Check10.p) annotation(
        Line(points = {{-38, 86}, {80, 86}, {80, -66}, {96, -66}, {96, -66}}, color = {0, 0, 127}));
      connect(pressure.y, eLC_Check9.p) annotation(
        Line(points = {{-38, 86}, {80, 86}, {80, -30}, {96, -30}, {96, -30}}, color = {0, 0, 127}));
      connect(pressure.y, eLC_Check6.p) annotation(
        Line(points = {{-38, 86}, {80, 86}, {80, 4}, {96, 4}, {96, 4}}, color = {0, 0, 127}));
      connect(pressure.y, eLC_Check7.p) annotation(
        Line(points = {{-38, 86}, {80, 86}, {80, 36}, {96, 36}, {96, 36}}, color = {0, 0, 127}));
      connect(pressure.y, eLC_Check8.p) annotation(
        Line(points = {{-38, 86}, {80, 86}, {80, 68}, {96, 68}, {96, 68}}, color = {0, 0, 127}));
      connect(pressure.y, eLC_Check5.p) annotation(
        Line(points = {{-38, 86}, {-20, 86}, {-20, -68}, {-4, -68}, {-4, -68}}, color = {0, 0, 127}));
      connect(pressure.y, eLC_Check4.p) annotation(
        Line(points = {{-38, 86}, {-20, 86}, {-20, -32}, {-4, -32}, {-4, -32}}, color = {0, 0, 127}));
      connect(pressure.y, eLC_Check3.p) annotation(
        Line(points = {{-38, 86}, {-20, 86}, {-20, 2}, {-4, 2}, {-4, 2}}, color = {0, 0, 127}));
      connect(pressure.y, eLC_Check2.p) annotation(
        Line(points = {{-38, 86}, {-20, 86}, {-20, 34}, {-4, 34}, {-4, 34}}, color = {0, 0, 127}));
      connect(pressure.y, eLC_Check1.p) annotation(
        Line(points = {{-38, 86}, {-20, 86}, {-20, 66}, {-4, 66}, {-4, 66}, {-4, 66}}, color = {0, 0, 127}));
      annotation(
        __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"));
    end ETC_Test2;

    model ETC_Test3
      IdealGasTest.EtaLamba.ELC_Check eLC_Check1(redeclare package Medium = Modelica.Media.IdealGases.SingleGases.Ar) annotation(
        Placement(visible = true, transformation(origin = {4, 62}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant pressure(k = 101325) annotation(
        Placement(visible = true, transformation(origin = {-50, 86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      IdealGasTest.EtaLamba.ELC_Check eLC_Check2(redeclare package Medium = Modelica.Media.IdealGases.SingleGases.CH4) annotation(
        Placement(visible = true, transformation(origin = {4, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      IdealGasTest.EtaLamba.ELC_Check eLC_Check3(redeclare package Medium = Modelica.Media.IdealGases.SingleGases.CO) annotation(
        Placement(visible = true, transformation(origin = {4, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      IdealGasTest.EtaLamba.ELC_Check eLC_Check4(redeclare package Medium = Modelica.Media.IdealGases.SingleGases.CO2) annotation(
        Placement(visible = true, transformation(origin = {4, -36}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      IdealGasTest.EtaLamba.ELC_Check eLC_Check5(redeclare package Medium = Modelica.Media.IdealGases.SingleGases.H2) annotation(
        Placement(visible = true, transformation(origin = {4, -72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      IdealGasTest.EtaLamba.ELC_Check eLC_Check6(redeclare package Medium = Modelica.Media.IdealGases.SingleGases.N2) annotation(
        Placement(visible = true, transformation(origin = {104, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      IdealGasTest.EtaLamba.ELC_Check eLC_Check7(redeclare package Medium = Modelica.Media.IdealGases.SingleGases.H2O) annotation(
        Placement(visible = true, transformation(origin = {104, 32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      IdealGasTest.EtaLamba.ELC_Check eLC_Check8(redeclare package Medium = Modelica.Media.IdealGases.SingleGases.O2) annotation(
        Placement(visible = true, transformation(origin = {104, 64}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      IdealGasTest.EtaLamba.ELC_Check eLC_Check9(redeclare package Medium = Modelica.Media.IdealGases.SingleGases.He) annotation(
        Placement(visible = true, transformation(origin = {104, -34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      IdealGasTest.EtaLamba.ELC_Check eLC_Check10(redeclare package Medium = Modelica.Media.IdealGases.SingleGases.C3H8) annotation(
        Placement(visible = true, transformation(origin = {104, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Ramp ramp1(duration = 1800, height = 1800, offset = 200, startTime = 0) annotation(
        Placement(visible = true, transformation(origin = {-50, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(ramp1.y, eLC_Check10.T) annotation(
        Line(points = {{-38, -80}, {66, -80}, {66, -74}, {96, -74}, {96, -74}}, color = {0, 0, 127}));
      connect(ramp1.y, eLC_Check9.T) annotation(
        Line(points = {{-38, -80}, {66, -80}, {66, -38}, {96, -38}, {96, -38}}, color = {0, 0, 127}));
      connect(ramp1.y, eLC_Check6.T) annotation(
        Line(points = {{-38, -80}, {66, -80}, {66, -4}, {96, -4}, {96, -4}}, color = {0, 0, 127}));
      connect(ramp1.y, eLC_Check7.T) annotation(
        Line(points = {{-38, -80}, {66, -80}, {66, 28}, {96, 28}, {96, 28}}, color = {0, 0, 127}));
      connect(ramp1.y, eLC_Check8.T) annotation(
        Line(points = {{-38, -80}, {66, -80}, {66, 60}, {96, 60}, {96, 60}}, color = {0, 0, 127}));
      connect(ramp1.y, eLC_Check5.T) annotation(
        Line(points = {{-38, -80}, {-30, -80}, {-30, -76}, {-4, -76}, {-4, -76}}, color = {0, 0, 127}));
      connect(ramp1.y, eLC_Check4.T) annotation(
        Line(points = {{-38, -80}, {-30, -80}, {-30, -40}, {-4, -40}, {-4, -40}}, color = {0, 0, 127}));
      connect(ramp1.y, eLC_Check3.T) annotation(
        Line(points = {{-38, -80}, {-30, -80}, {-30, -6}, {-4, -6}, {-4, -6}}, color = {0, 0, 127}));
      connect(ramp1.y, eLC_Check2.T) annotation(
        Line(points = {{-38, -80}, {-30, -80}, {-30, 26}, {-4, 26}, {-4, 26}}, color = {0, 0, 127}));
      connect(ramp1.y, eLC_Check1.T) annotation(
        Line(points = {{-38, -80}, {-30, -80}, {-30, 58}, {-4, 58}, {-4, 58}}, color = {0, 0, 127}));
      connect(pressure.y, eLC_Check10.p) annotation(
        Line(points = {{-38, 86}, {80, 86}, {80, -66}, {96, -66}, {96, -66}}, color = {0, 0, 127}));
      connect(pressure.y, eLC_Check9.p) annotation(
        Line(points = {{-38, 86}, {80, 86}, {80, -30}, {96, -30}, {96, -30}}, color = {0, 0, 127}));
      connect(pressure.y, eLC_Check6.p) annotation(
        Line(points = {{-38, 86}, {80, 86}, {80, 4}, {96, 4}, {96, 4}}, color = {0, 0, 127}));
      connect(pressure.y, eLC_Check7.p) annotation(
        Line(points = {{-38, 86}, {80, 86}, {80, 36}, {96, 36}, {96, 36}}, color = {0, 0, 127}));
      connect(pressure.y, eLC_Check8.p) annotation(
        Line(points = {{-38, 86}, {80, 86}, {80, 68}, {96, 68}, {96, 68}}, color = {0, 0, 127}));
      connect(pressure.y, eLC_Check5.p) annotation(
        Line(points = {{-38, 86}, {-20, 86}, {-20, -68}, {-4, -68}, {-4, -68}}, color = {0, 0, 127}));
      connect(pressure.y, eLC_Check4.p) annotation(
        Line(points = {{-38, 86}, {-20, 86}, {-20, -32}, {-4, -32}, {-4, -32}}, color = {0, 0, 127}));
      connect(pressure.y, eLC_Check3.p) annotation(
        Line(points = {{-38, 86}, {-20, 86}, {-20, 2}, {-4, 2}, {-4, 2}}, color = {0, 0, 127}));
      connect(pressure.y, eLC_Check2.p) annotation(
        Line(points = {{-38, 86}, {-20, 86}, {-20, 34}, {-4, 34}, {-4, 34}}, color = {0, 0, 127}));
      connect(pressure.y, eLC_Check1.p) annotation(
        Line(points = {{-38, 86}, {-20, 86}, {-20, 66}, {-4, 66}, {-4, 66}, {-4, 66}}, color = {0, 0, 127}));
      annotation(
        __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"),
        experiment(StartTime = 0, StopTime = 1800, Tolerance = 1e-6, Interval = 10));
    end ETC_Test3;










  end EtaLamba;

  package DryAirNasaTest
    model CheckEL
      replaceable package Medium = Modelica.Media.Air.DryAirNasa;
      Medium.ThermodynamicState state;
      parameter Medium.AbsolutePressure p = 101325.0;
      Modelica.Blocks.Interfaces.RealInput T annotation(
        Placement(visible = true, transformation(origin = {-74, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-74, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput eta annotation(
        Placement(visible = true, transformation(origin = {86, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {84, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput lambda annotation(
        Placement(visible = true, transformation(origin = {84, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {86, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      state = Medium.setState_pT(p, T);
      eta = Medium.dynamicViscosity(state);
      lambda = Medium.thermalConductivity(state);
      annotation(
        Icon(graphics = {Text(origin = {-49, 9}, extent = {{1, 5}, {21, -33}}, textString = "T"), Text(origin = {34, 67}, extent = {{-20, 11}, {32, -23}}, textString = "eta"), Text(origin = {39, -76}, extent = {{-77, 54}, {29, -22}}, textString = "lambda"), Rectangle(origin = {0, 1}, extent = {{-80, 79}, {80, -81}})}));
    end CheckEL;

    model CheckELTest
    IdealGasTest.DryAirNasaTest.CheckEL checkEL1 annotation(
        Placement(visible = true, transformation(origin = {22, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Ramp temperature(duration = 1700, height = 1700, offset = 200, startTime = 0)  annotation(
        Placement(visible = true, transformation(origin = {-34, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(temperature.y, checkEL1.T) annotation(
        Line(points = {{-22, 0}, {14, 0}, {14, 0}, {14, 0}}, color = {0, 0, 127}));
    annotation(
        experiment(StartTime = 0, StopTime = 1700, Tolerance = 1e-6, Interval = 50),
        __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"));end CheckELTest;




  end DryAirNasaTest;
  annotation(
    uses(Modelica(version = "3.2.2")));
end IdealGasTest;