'use strict';

angular.module('clotho.trails')
	.controller('TestTrailSplashCtrl', function ($scope, $location, Clotho) {
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

		 trail.materials = [{
		 "name" : "Course Slides",
		 "icon" : "syllabus",
		 "type" : "Syllabus",
		 "path" : "partials/trails/clotho-intro/materials/ClothoIntroSyllabus.pdf"
		 }]

		//to get all the trail IDs above

		 var trails = [];
		 angular.forEach($scope.topics, function (topic) {
		  angular.forEach(topic.trails, function (trail) {
		    trails.push(trail.id)
		  });
		 });


		as of 5/22/14, but probably not worth re-importing after all the customizations made

		 var list = ["PL2aPXzks-TgOacBZ1TGGwANuld5yGmbkt", "PL2aPXzks-TgOj55WalLbDUTfmuhZvYVmW", "PL2aPXzks-TgNK702CfJaz-gTfRtv34sNt", "PL2aPXzks-TgMRN9PNtntlYoMuV5Xb8iu2", "PL2aPXzks-TgMvy4PG5wwbrw9VKBuowOR8", "PL2aPXzks-TgPIVgTQU-ajK2eCP9deCIl2", "PL2aPXzks-TgPX7HU4pogNUSBN9OK2PhBK", "PL2aPXzks-TgNkKi2lUWL65MpPsRBCb5q5", "PL2aPXzks-TgPagNtoxbU5Kqt1aZCm1ORd", "PL2aPXzks-TgPiVtQqqthOAE60E7iHz-gg", "PL2aPXzks-TgOOPfxneYJ3pBwTlqzv9Ezn", "PL2aPXzks-TgN9w3Pm49zEQJfnNDV8T1KZ", "PL2aPXzks-TgNwbIY4SU7Wd1RObYSllKqv", "PL2aPXzks-TgMgibxKB80naduiRZnqcZPR", "PL2aPXzks-TgOtngzKdocz4P50TCncRyPy", "PL2aPXzks-TgNn3gT7KNn3m2QN3x_m4NCW", "PL2aPXzks-TgN_Ztp-bRkzbtbi5Ncsh71t", "PL2aPXzks-TgPrTigCFi0P8hfWz-CEO5qL", "PL2aPXzks-TgPYH_y0UNKluRkrRHbcQLuX", "PL2aPXzks-TgNlrgM8BK3jD_aiyjoXOsRI", "PL2aPXzks-TgPMLXx10saMCO7oKp7IMKDn", "PL2aPXzks-TgMJr2fAKj5ixn9RXypJ_zaZ", "PL2aPXzks-TgPhvEG_0c44SFt03b5sUGov", "PL2aPXzks-TgNe-LNrcRmB-iHvFFZokK7R", "PL2aPXzks-TgNmfjrolgYWGblQDNzo11Yq", "PL2aPXzks-TgMbE-b15ezcHQmGVTfMuHcZ", "PL2aPXzks-TgNu5fDsnbCOdQqdSqU6m1FF", "PL2aPXzks-TgPkmcVYYZACTbbq4bx_ZAoI", "PL2aPXzks-TgNP1FuAv0v-fS6UhCjy1waN", "PL2aPXzks-TgP45-xudcn4XXiVZfUOYMht", "PL2aPXzks-TgMz06rZlbCRKzYyVkXyIwb9", "PL2aPXzks-TgPQbhJWXI3fSsV9FZ0QLAXa", "PL2aPXzks-TgOdmDmfWSneeIdp_4jF9lI-", "PL2aPXzks-TgP0KMKfq4-fQPKClxvkYwWd", "PL2aPXzks-TgO0k9PhT__NSh2x6HNimaOy", "PL2aPXzks-TgOYRv5T9yRzwJE1QWEiY2R8", "PL2aPXzks-TgMo_x6XhYmfoZDezm4m-DcL"]


		$scope.createAll = function () {
			angular.forEach(list, function (playlistid) {
				Youtube.playlistToTrail(playlistid).then(function (result) {
					Clotho.create(result);
				});
			});
		};
		*/

		$scope.startTrail = function (id) {
			console.log(id);
			Clotho.startTrail(id);
		};

		$scope.startTrailPage = function (id, pos) {
			$location.search('position', pos);
			Clotho.startTrail(id);
		};

		$scope.highlight = function (trail, evt) {
			//highlight is trail in object above, selected is the actual trail
			$scope.highlighted = trail;
			$scope.loading = true;
			Clotho.get(trail.id, {mute : true}).then(function (result) {
				$scope.loading = false;
				$scope.selected = result;
			});
		}
	});
