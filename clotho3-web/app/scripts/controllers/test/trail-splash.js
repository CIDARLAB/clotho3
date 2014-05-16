'use strict';

angular.module('clotho.trails')
	.controller('TestTrailSplashCtrl', function ($scope, Youtube) {
		$scope.topics = [
			{
				"title": "Introducing Clotho",
				"trails": [
					{ title: "Introduction to Clotho", "playlist": "" }
				]
			},
			{
				"title": "Basic Molecular Biology",
				"trails": [
					{ title: "Introduction to Synthetic Biology", "playlist": "PL2aPXzks-TgPYH_y0UNKluRkrRHbcQLuX" },
					{ title: "Organic Chemistry and Molecular Biology", "playlist": "PL2aPXzks-TgMbE-b15ezcHQmGVTfMuHcZ" }
				]
			},
			{
				"title": "Fabrication",
				"trails": [
					{ title: "DNA Manipulation Enzymes", "playlist": "PL2aPXzks-TgNkKi2lUWL65MpPsRBCb5q5" },
					{ title: "DNA Fabrication", "playlist": "PL2aPXzks-TgPX7HU4pogNUSBN9OK2PhBK" },
					{ title: "Genome Manipulation", "playlist": "PL2aPXzks-TgOtngzKdocz4P50TCncRyPy" },
					{ title: "Chassis and Strains", "playlist": "PL2aPXzks-TgMvy4PG5wwbrw9VKBuowOR8" },
					{ title: "Parts", "playlist": "PL2aPXzks-TgNP1FuAv0v-fS6UhCjy1waN" },
					{ title: "Human Practices", "playlist": "PL2aPXzks-TgN_Ztp-bRkzbtbi5Ncsh71t" }
				]
			},
			{
				"title": "Analysis",
				"trails": [
					{ title: "Whole-Cell Analysis", "playlist": "PL2aPXzks-TgMo_x6XhYmfoZDezm4m-DcL" },
					{ title: "In vitro Analysis", "playlist": "PL2aPXzks-TgNlrgM8BK3jD_aiyjoXOsRI" },
					{ title: "Directed Evolution", "playlist": "PL2aPXzks-TgOOPfxneYJ3pBwTlqzv9Ezn" },
					{ title: "Combinatorial Libraries", "playlist": "PL2aPXzks-TgPIVgTQU-ajK2eCP9deCIl2" },
					{ title: "Automation Tools", "playlist": "PL2aPXzks-TgOj55WalLbDUTfmuhZvYVmW" }
				]
			},
			{
				"title": "Design",
				"trails": [
					{ title: "Basics of OOP and Organic Chemistry", "playlist": "PL2aPXzks-TgOv8a6ycw13PCY0UcsTxNkA" },
					{ title: "Introduction to Biosynthesis", "playlist": "PL2aPXzks-TgPrTigCFi0P8hfWz-CEO5qL" },
					{ title: "Literature Search and Pathway Discovery", "playlist": "PL2aPXzks-TgMJr2fAKj5ixn9RXypJ_zaZ" },
					{ title: "Metabolic Optimization", "playlist": "PL2aPXzks-TgPhvEG_0c44SFt03b5sUGov" },
					{ title: "Monofunctional Enzyme Biosynthesis Examples ", "playlist": "PL2aPXzks-TgNmfjrolgYWGblQDNzo11Yq" },
					{ title: "Polymeric Metabolites", "playlist": "PL2aPXzks-TgP45-xudcn4XXiVZfUOYMht" },
					{ title: "Modular PKS and NRP Engineering", "playlist": "PL2aPXzks-TgNe-LNrcRmB-iHvFFZokK7R" },
					{ title: "Processes and Localization", "playlist": "PL2aPXzks-TgMz06rZlbCRKzYyVkXyIwb9" },
					{ title: "Transcription Models", "playlist": "PL2aPXzks-TgO0k9PhT__NSh2x6HNimaOy" },
					{ title: "Translation Models", "playlist": "PL2aPXzks-TgOYRv5T9yRzwJE1QWEiY2R8" },
					{ title: "Biochemical Circuits", "playlist": "PL2aPXzks-TgNK702CfJaz-gTfRtv34sNt" },
					{ title: "Stochastic modeling", "playlist": "PL2aPXzks-TgP0KMKfq4-fQPKClxvkYwWd" },
					{ title: "Biological Relationships", "playlist": "PL2aPXzks-TgMRN9PNtntlYoMuV5Xb8iu2" },
					{ title: "Engineering microbial interactions", "playlist": "PL2aPXzks-TgN9w3Pm49zEQJfnNDV8T1KZ" },
					{ title: "Refactoring", "playlist": "PL2aPXzks-TgPQbhJWXI3fSsV9FZ0QLAXa" },
					{ title: "Advanced microbial processes", "playlist": "PL2aPXzks-TgOacBZ1TGGwANuld5yGmbkt" },
					{ title: "Eukaryotic Devices", "playlist": "PL2aPXzks-TgNwbIY4SU7Wd1RObYSllKqv" },
					{ title: "Developmental Devices", "playlist": "PL2aPXzks-TgPagNtoxbU5Kqt1aZCm1ORd" },
					{ title: "Part Engineering", "playlist": "PL2aPXzks-TgPkmcVYYZACTbbq4bx_ZAoI" }
				]
			},
			{
				"title": "Wetlab: Cloning Experiment",
				"trails": [
					{ title: "Overview of Basic Cloning", "playlist": "PL2aPXzks-TgNu5fDsnbCOdQqdSqU6m1FF" },
					{ title: "Running PCR", "playlist": "PL2aPXzks-TgOdmDmfWSneeIdp_4jF9lI-" },
					{ title: "Digesting DNAs", "playlist": "PL2aPXzks-TgPiVtQqqthOAE60E7iHz-gg" },
					{ title: "Gel Purification", "playlist": "PL2aPXzks-TgMgibxKB80naduiRZnqcZPR" },
					{ title: "Ligation and Transformation", "playlist": "PL2aPXzks-TgPMLXx10saMCO7oKp7IMKDn" },
					{ title: "Growing up colonies and prepping DNA", "playlist": "PL2aPXzks-TgNn3gT7KNn3m2QN3x_m4NCW" }
				]
			},
			{
				"title": "Advanced Tools",
				"trails": [

				]
			}
		];

		$scope.startTrail = function (trail) {
			//todo
		};

		$scope.highlight = function (trail) {
			$scope.highlighted = trail;
			//note - in interim, let's just fetch the playlist and construct lazily
			Youtube.playlistToTrail(trail.playlist).then(function (compiled) {
				$scope.selected = compiled;
			});
		}
	});
