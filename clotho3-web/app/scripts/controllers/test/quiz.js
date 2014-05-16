'use strict';

angular.module('clotho.webapp')
	.controller('TestQuizCtrl', function ($scope) {

		// todo - experiment with lazily retrieved quizzes

		// todo - example self-contained grading function

		// fixme - change dynamic dictionary syntax - not sloppy
		// todo - refactor from names to IDs for get() and run()

		$scope.trueFalse = {
			question: {
				type: "truefalse",
				title: "True False",
				question: "Clotho is great?"
			},
			grade: {
				answer: {
					type: "boolean",
					value: true
				}
			}
		};

		$scope.number = {
			question: {
				type: "number",
				title: "Number with tolerance (0.01)",
				question: "What is Pi?"
			},
			grade: {
				answer: {
					type: "number",
					tolerance: 0.01,
					value: 3.1415926
				}
			}
		};

		$scope.staticMC = {
			question: {
				type: "mc",
				title: "Static Multiple Choice",
				question: "Find the reverse complement of the sequence ACCGGGTTTT",
				options: ["accgggtttt", "TGGCCCAAAA", "TTTTGGGCCA", "AAAACCCGGT"],
				hint: "Example hint!"
			},
			options: {
				showAnswer: true,
				allowMultiple: true,
				allowRetry: true,
				randomization: false
			},
			grade: {
				answer: {
					type: "string",
					value: "AAAACCCGGT"
				}
			}
		};

		$scope.staticTemplating = {
			question: {
				type: "mc",
				title: "Static Templating",
				question: "Find the reverse complement of the sequence {{mySeq}}",
				options: ["{{value1}}", "{{value2}}", "{{value3}}", "{{value4}}"]
			},
			dictionary: {
				static: {
					mySeq: "ACCGGGTTTT",
					value1: "accgggtttt",
					value2: "TGGCCCAAAA",
					value3: "TTTTGGGCCA",
					value4: "AAAACCCGGT"
				}
			},
			grade: {
				answer: {
					type: "string",
					value: "AAAACCCGGT"
				}
			}
		};

		$scope.staticFeedback = {
			question: {
				type: "mc",
				title: "Static Feedback",
				question: "Find the reverse complement of the sequence ACCGGGTTTT",
				options: ["accgggtttt", "TGGCCCAAAA", "TTTTGGGCCA", "AAAACCCGGT"]
			},
			options: {
				allowRetry: true
			},
			grade: {
				answer: {
					type: "string",
					value: "AAAACCCGGT"
				}
			},
			feedback: {
				default: 'To reverse complement a sequence, simply reverse the order, and complement each base (order in which you do this does not matter)',
				static: {
					"accgggtttt": "Nope, that's the lowercase",
					"TGGCCCAAAA": "Nope, that's the complement",
					"TTTTGGGCCA": "Nope, that's the reverse"
				}
			}
		};

		$scope.staticFillin = {
			question: {
				type: "fillin",
				title: "Static Fillin",
				question: "Find the reverse complement of the sequence AAACCGT"
			},
			options: {
				showAnswer: true
			},
			grade: {
				answer: {
					type: "string",
					value: "ACGGTTT"
				}
			}
		};

		$scope.dynamicTemplating = {
			question: {
				type: "mc",
				title: "Dynamic Templating",
				question: "Find the reverse complement of the sequence {{mySeq}}",
				options: ["{{value1}}", "{{value2}}", "{{value3}}", "{{mySeq}}"],
				hint: "When you retry, all templated values are recalculated. They are then passed to the server for grading"
			},
			options: {
				allowRetry: true
			},
			dictionary: {
				dynamic: [
					{
						mySeq: {
							id : "org.clothocad.test.randomSequence",
							args : [16]
						}
					},
					{
						value1: {
							id : "org.clothocad.test.revcomp",
							args : ["{{mySeq}}"]
						}
					},
					{
						value2: {
							id : "org.clothocad.test.reverse",
							args : ["{{mySeq}}"]
						}
					},
					{
						value3: {
							id : "org.clothocad.test.complement",
							args : ["{{mySeq}}"]
						}
					}
				]
			},
			grade: {
				answer: {
					type: "function",
					value: "org.clothocad.test.revcomp"
				},
				args: ['{{mySeq}}']
			}
		};

		$scope.functionGrading = {
			question: {
				type: "fillin",
				title: "Own grading function",
				question: "Copy in this sequence: {{mySeq}}",
				hint: "This question has its own grading function, for more complicated logic. Arguments can be specified by author. This one just makes sure args + input are equal."
			},
			options: {
				allowRetry: true
			},
			dictionary: {
				dynamic: [
					{
						mySeq: {
							id : "org.clothocad.test.randomSequence",
							args : [40]
						}
					}
				]
			},
			grade: {
				function: "org.clothocad.test.customGradeFunction",
				args: ['{{mySeq}}']
			}
		};

	});
