    Info<< "Reading field T\n" << endl;
    
    volScalarField Rate
    (
        IOobject
        (
            "Rate",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField B
    (
        IOobject
        (
            "B",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    volScalarField A
    (
        IOobject
        (
            "A",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    volScalarField C
    (
        IOobject
        (
            "C",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );




    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );


    Info<< "Reading diffusivity DA\n" << endl;

    dimensionedScalar DA
    (
        transportProperties.lookup("DA")
    );
    
        Info<< "Reading diffusivity DB\n" << endl;

    dimensionedScalar DB
    (
        transportProperties.lookup("DB")
    );
    
    dimensionedScalar k
    (
        transportProperties.lookup("k")
    );

 
    dimensionedScalar gamma
    (
        transportProperties.lookup("gamma")
    );



