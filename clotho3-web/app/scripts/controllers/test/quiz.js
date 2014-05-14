'use strict';

angular.module('clotho.webapp')
  .controller('TestQuizCtrl', function ($scope) {

		//todo - refactor from names to IDs for get() and run()

		$scope.staticMC = {
			question : {
				type : "mc",
				question : "Find the reverse complement of the sequence ACCGGGTTTT",
				options : ["accgggtttt", "TGGCCCAAAA", "TTTTGGGCCA", "AAAACCCGGT"]
			},
			options : {
				checkAnswer : true,
				showAnswer : true,
				allowMultiple : true,
				allowRetry : true,
				randomization : false
			},
			grade: {
				answer : "AAAACCCGGT"
			}
		},

		$scope.staticTemplating = {
			question : {
				type : "mc",
				question : "Find the reverse complement of the sequence {{mySeq}}",
				options : ["{{value1}}", "{{value2}}", "{{value2}}", "{{value2}}"]
			},
			options : {
				multipleAttempts : true,
				retry : true
			},
			dictionary: {
				static : [
					{ mySeq : "ACCGGGTTTT" },
					{ value1 : "accgggtttt" },
					{ value2 : "TGGCCCAAAA" },
					{ value3 : "TTTTGGGCCA" },
					{ value4 : "AAAACCCGGT" }
				]
			},
			grade: {
				answer : "AAAACCCGGT"
			}
		};

		$scope.staticFeedback = {
			question : {
				type : "mc",
				question : "Find the reverse complement of the sequence ACCGGGTTTT",
				options : ["accgggtttt", "TGGCCCAAAA", "TTTTGGGCCA", "AAAACCCGGT"]
			},
			grade: {
				answer : "AAAACCCGGT"
			},
			feedback : {
				default : 'Feedback shown when other keys are not met',
				static : [
					{ "accgggtttt" : "Nope, that's the lowercase" },
					{ "TGGCCCAAAA" : "Nope, that's the complement" },
					{ "TTTTGGGCCA" : "Nope, that's the reverse" }
				]
			}
		};

		$scope.staticFillin = {
			question : {
				type : "fillin",
				question : "Find the reverse complement of the sequence AAACCGT"
			},
			grade: {
				answer : {
					type : "string",
					value : "ACGGTTT"
				}
			}
		};

		$scope.dynamicTemplating = {
			question : {
				type : "mc",
				question : "Find the reverse complement of the sequence {{mySeq}}",
				options : ["{{value1}}", "{{value2}}", "{{value3}}", "{{value4}}"]
			},
			dictionary: {
				dynamic : [
					{mySeq  : "clotho.run('randomSequence', [16])" },
					{value1 : "clotho.run('DNA.revcomp', ['{{mySeq}}'])" },
					{value2 : "clotho.run('DNA.reverse', ['{{mySeq}}'])" },
					{value3 : "clotho.run('DNA.complement', ['{{mySeq}}'])" },
					{value4 : "clotho.run('reverse', ['{{mySeq}}'])" }
				]
			},
			grade: {
				function : "DNA.revcomp",
				args : ['{{mySeq}}']
			}
		}

  });
