//- Create Objects for calculation
IdealReactorProperties properties("");

Thermo thermo(properties.thermo());

Transport transport(properties.transport(), thermo);

Chemistry chemistry(properties.chemistry(), thermo);

//- Interprete data and store for analysis in files
if (properties.interprete())
{
    Interpreter interpreter;

    interpreter.summary(transport, thermo, chemistry);

    Footer(startTime);
    return 0;
}
