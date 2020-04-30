//- Create Objects for calculation
IdealReactorProperties properties("");

Thermo thermo(properties.thermo());

Transport transport(properties.transport(), thermo);

Chemistry chemistry(properties.chemistry(), thermo);

//- Insert species word list to the Transport object
transport.insertChemistrySpecies(chemistry.species());

//- Interprete data and store for analysis in files
if (properties.interprete())
{
    Interpreter interpreter;

    interpreter.summary(transport, thermo, chemistry);

    Footer(startTime);
    return 0;
}

//- Create Time object
Time time(properties.dict());
