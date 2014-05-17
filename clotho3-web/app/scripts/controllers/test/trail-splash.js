'use strict';

angular.module('clotho.trails')
	.controller('TestTrailSplashCtrl', function ($scope, Clotho, Youtube) {
		$scope.topics = [
			{
				"title": "Introducing Clotho",
				"trails": [
					{
						"title": "Introduction to Clotho",
						"id": "org.clothocad.trails.LearningClotho"
					}
				]
			},
			{
				"title": "Basic Molecular Biology",
				"trails": [
					{
						"title": "Introduction to Synthetic Biology",
						"id": "org.clothocad.trails.youtube.IntroductiontoSyntheticBiology"
					},
					{
						"title": "Organic Chemistry and Molecular Biology",
						"id": "org.clothocad.trails.youtube.OrganicChemistryandMolecularBiology"
					}
				]
			},
			{
				"title": "Fabrication",
				"trails": [
					{
						"title": "DNA Manipulation Enzymes",
						"id": "org.clothocad.trails.youtube.DNAManipulationEnzymes"
					},
					{
						"title": "DNA Fabrication",
						"id": "org.clothocad.trails.youtube.DNAFabrication"
					},
					{
						"title": "Genome Manipulation",
						"id": "org.clothocad.trails.youtube.GenomeManipulation"
					},
					{
						"title": "Chassis and Strains",
						"id": "org.clothocad.trails.youtube.ChassisandStrains"
					},
					{
						"title": "Parts",
						"id": "org.clothocad.trails.youtube.Parts"
					},
					{
						"title": "Human Practices",
						"id": "org.clothocad.trails.youtube.HumanPractices"
					}
				]
			},
			{
				"title": "Analysis",
				"trails": [
					{
						"title": "Whole-Cell Analysis",
						"id": "org.clothocad.trails.youtube.Whole-CellAnalysis"
					},
					{
						"title": "In vitro Analysis",
						"id": "org.clothocad.trails.youtube.InvitroAnalysis"
					},
					{
						"title": "Directed Evolution",
						"id": "org.clothocad.trails.youtube.DirectedEvolution"
					},
					{
						"title": "Combinatorial Libraries",
						"id": "org.clothocad.trails.youtube.CombinatorialLibraries"
					},
					{
						"title": "Automation Tools",
						"id": "org.clothocad.trails.youtube.AutomationTools"
					}
				]
			},
			{
				"title": "Design",
				"trails": [
					{
						"title": "Introduction to Biosynthesis",
						"id": "org.clothocad.trails.youtube.IntroductiontoBiosynthesis"
					},
					{
						"title": "Literature Search and Pathway Discovery",
						"id": "org.clothocad.trails.youtube.LiteratureSearchandPathwayDiscovery"
					},
					{
						"title": "Metabolic Optimization",
						"id": "org.clothocad.trails.youtube.MetabolicOptimization"
					},
					{
						"title": "Monofunctional Enzyme Biosynthesis Examples ",
						"id": "org.clothocad.trails.youtube.MonofunctionalEnzymeBiosynthesisExamples"
					},
					{
						"title": "Polymeric Metabolites",
						"id": "org.clothocad.trails.youtube.PolymericMetabolites"
					},
					{
						"title": "Modular PKS and NRP Engineering",
						"id": "org.clothocad.trails.youtube.ModularPKSandNRPEngineering"
					},
					{
						"title": "Processes and Localization",
						"id": "org.clothocad.trails.youtube.ProcessesandLocalization"
					},
					{
						"title": "Transcription Models",
						"id": "org.clothocad.trails.youtube.TranscriptionModels"
					},
					{
						"title": "Translation Models",
						"id": "org.clothocad.trails.youtube.TranslationModels"
					},
					{
						"title": "Biochemical Circuits",
						"id": "org.clothocad.trails.youtube.BiochemicalCircuits"
					},
					{
						"title": "Stochastic modeling",
						"id": "org.clothocad.trails.youtube.Stochasticmodeling"
					},
					{
						"title": "Biological Relationships",
						"id": "org.clothocad.trails.youtube.BiologicalRelationships"
					},
					{
						"title": "Engineering microbial interactions",
						"id": "org.clothocad.trails.youtube.Engineeringmicrobialinteractions"
					},
					{
						"title": "Refactoring",
						"id": "org.clothocad.trails.youtube.Refactoring"
					},
					{
						"title": "Advanced microbial processes",
						"id": "org.clothocad.trails.youtube.Advancedmicrobialprocesses"
					},
					{
						"title": "Eukaryotic Devices",
						"id": "org.clothocad.trails.youtube.EukaryoticDevices"
					},
					{
						"title": "Developmental Devices",
						"id": "org.clothocad.trails.youtube.DevelopmentalDevices"
					},
					{
						"title": "Part Engineering",
						"id": "org.clothocad.trails.youtube.PartEngineering"
					}
				]
			},
			{
				"title": "Wetlab: Cloning Experiment",
				"trails": [
					{
						"title": "Overview of Basic Cloning",
						"id": "org.clothocad.trails.youtube.OverviewofBasicCloning"
					},
					{
						"title": "Running PCR",
						"id": "org.clothocad.trails.youtube.RunningPCR"
					},
					{
						"title": "Digesting DNAs with restriction endonucleases",
						"id": "org.clothocad.trails.youtube.DigestingDNAswithrestrictionendonucleases"
					},
					{
						"title": "Gel Purification",
						"id": "org.clothocad.trails.youtube.GelPurification"
					},
					{
						"title": "Ligation and Transformation",
						"id": "org.clothocad.trails.youtube.LigationandTransformation"
					},
					{
						"title": "Growing up colonies and then prepping DNA",
						"id": "org.clothocad.trails.youtube.GrowingupcoloniesandthenpreppingDNA"
					}
				]
			},
			{
				"title": "Advanced Tools",
				"trails": []
			}
		];

		/*
		$scope.createAll = function () {
			angular.forEach($scope.topics, function (topic) {
				angular.forEach(topic.trails, function (trail) {
					Youtube.playlistToTrail(trail.playlist).then(function (result) {
						Clotho.create(result);
					});
				});
			});
		};
		*/

		$scope.startTrail = function (trail) {
			//todo
		};

		$scope.highlight = function (trail, evt) {
			$scope.highlighted = trail;
			$scope.loading = true;
			//note - in interim, let's just fetch the playlist and construct lazily
			//Youtube.playlistToTrail(trail.playlist).then(function (compiled) {
			Clotho.get(trail.id).then(function (result) {
				//cna't use translate because not consistent across browsers (vendor prefixing)
				$scope.highlightStyle = {
					top: evt.target.offsetTop
				};
				$scope.loading = false;
				$scope.selected = result;
			});
		}
	});
