'use strict';

angular.module('clotho.webapp')
  .controller('TestConstructionCtrl', function ($scope) {

		$scope.parsed = {
			"name": "BioBricking of Kanamycin Gene",
			"schema": "org.clothocad.model.constructionFile",
			"description": "BioBricks is a standard assembly using a defined set of enzymes, enabling easy, directional addition of parts",
			"steps": [
				{
					"reaction": "pcr",
					"input": [
						"pSB1AK3-b0015",
						[
							"ca1067F",
							"ca1067R"
						]
					],
					"output": "pcrpdt"
				},
				{
					"reaction": "digest",
					"input": [
						"pcrpdt",
						"No enzymes specified"
					],
					"sizes": "1038+10+6",
					"choice": "L",
					"output": "pcrdig"
				},
				{
					"reaction": "digest",
					"input": [
						"pSB1A2-I13521",
						[
							"EcoRI",
							"SpeI"
						]
					],
					"sizes": "2062+946",
					"choice": "L",
					"output": "vectdig"
				},
				{
					"reaction": "ligate",
					"input": [
						[
							"pcrdig",
							"vectdig"
						]
					],
					"output": "pSB1A2-Bca9128"
				}
			],
			"dictionary": [
				{
					"key": "ca1067f",
					"description": "Biobricking of KanR of pSB1AK3",
					"value": "ccagtGAATTCgtccTCTAGAgagctgatccttcaactc"
				},
				{
					"key": "ca1067r",
					"description": "Biobricking of KanR of pSB1AK3",
					"value": "gcagtACTAGTtccgtcaagtcagcgtaatg"
				}
			]
		};


	});
