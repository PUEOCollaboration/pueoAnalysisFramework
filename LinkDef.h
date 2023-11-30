#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace pueo; 

#pragma link C++ class pueo::EventReconstructor+;
#pragma link C++ class pueo::EventSummary+;
#pragma link C++ class pueo::EventSummary::PointingHypothesis+;
#pragma link C++ class pueo::EventSummary::WaveformInfo+;
#pragma link C++ class pueo::EventSummary::ChannelInfo+;
#pragma link C++ class pueo::EventSummary::EventFlags+;
#pragma link C++ class pueo::EventSummary::SourceHypothesis+;
#pragma link C++ class pueo::EventSummary::MCTruth+;
#pragma link C++ class pueo::EventSummary::PayloadLocation+;
#pragma link C++ class pueo::FilteredEvent+;
#pragma link C++ class pueo::FilterStrategy+;
#pragma link C++ class pueo::AnalysisWaveform+;
#pragma link C++ class pueo::AnalysisWaveform::PowerCalculationOptions;
#pragma link C++ class pueo::PrettyAnalysisWaveform+;
#pragma link C++ class pueo::TGraphAligned+;

#pragma link C++ class pueo::TemplateSummary+;
#pragma link C++ class pueo::TemplateSummary::SingleTemplateResult+;
#pragma link C++ class pueo::TemplateMachine+;

#pragma link C++ class pueo::CorrelationSummary+;
#pragma link C++ class pueo::NoiseSummary+;
#pragma link C++ class pueo::NoiseMonitor+;

#pragma link C++ class pueo::FilterOperation+;
#pragma link C++ class pueo::UniformFilterOperation+;
#pragma link C++ class pueo::ConditionalFilterOperation+;
#pragma link C++ class pueo::SimplePassBandFilter+;
#pragma link C++ class pueo::SimpleNotchFilter+;
#pragma link C++ class pueo::HybridFilter+;
#pragma link C++ class pueo::SumDifferenceFilter+;
#pragma link C++ class pueo::DigitalFilterOperation+;
#pragma link C++ class pueo::GeometricFilter+;
#pragma link C++ class pueo::GaussianTaper; 
#pragma link C++ class pueo::DiodeFilter; 
#pragma link C++ class pueo::DeglitchFilter+;
#pragma link C++ class pueo::InterpolatedNotchFilter; 

#pragma link C++ namespace pueo::Response+;
#pragma link C++ namespace pueo::impulsivity+;
#pragma link C++ namespace pueo::bandwidth+;

#pragma link C++ namespace pueo::polarimetry+;
#pragma link C++ class pueo::polarimetry::StokesAnalysis;

#pragma link C++ class pueo::DeconvolutionMethod+;
#pragma link C++ class pueo::NaiveDeconvolution+;
#pragma link C++ class pueo::BandLimitedDeconvolution+;
#pragma link C++ class pueo::CLEANDeconvolution+;
#pragma link C++ class pueo::AbstractResponse+;
#pragma link C++ class pueo::Response+;
#pragma link C++ class pueo::ResponseManager+;
#pragma link C++ class pueo::CompositeResponse+;

#pragma link C++ class pueo::DeconvolveFilter+;
#pragma link C++ class pueo::WienerDeconvolution+; 
#pragma link C++ class pueo::AllPassDeconvolution+; 
#pragma link C++ class pueo::ImpulseResponseXCorr+; 
#pragma link C++ class pueo::CLEAN+; 
#pragma link C++ class pueo::DeconvolutionMethod+;

#pragma link C++ class pueo::EventFaker+; 
#pragma link C++ class pueo::SensitivityCalculator+; 
#pragma link C++ class pueo::FreqDomainFunction+; 
#pragma link C++ class pueo::PayloadParameters+; 
#pragma link C++ class pueo::UsefulAttitude+; 


#endif


