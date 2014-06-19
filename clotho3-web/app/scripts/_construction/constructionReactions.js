angular.module('clotho.construction')
/**
 * @ngdoc value
 * @name ConstructionReactions
 * @description map of reactions to Clotho Function IDs and other info about them. Add to this list as new ones are supported
 */
	.value('ConstructionReactions', {
		"pcr": {
			reactionId: "clotho.functions.dna.pcr",
			readable: "PCR",
			template_step : "views/_construction/step/pcr.html",
			template_wetlabL : "PCR {{step.input[1]}}/{{step.input[2]}} on {{step.input[0]}}",
			template_wetlabR : "({{step.output}})"
		},
		"ligate": {
			reactionId: "clotho.functions.dna.ligate",
			readable: "Ligate",
			template_step : "views/_construction/step/ligate.html",
			template_wetlabL : "Ligate {{step.input[0] | joinArray:'/'}}",
			template_wetlabR : "({{step.output}})"
		},
		"digest": {
			reactionId: "clotho.functions.dna.digest",
			readable: "Digest",
			template_step : "views/_construction/step/digest.html",
			template_wetlabL : "Digest {{step.input[0]}}",
			template_wetlabR : "({{step.input[1] | joinArray:'/'}}, {{step.output}})"
		},
		"gelpurify": {
			reactionId: "clotho.functions.dna.gelPurify",
			readable: "Gel Purify",
			template_step : "views/_construction/step/gelpurify.html",
			template_wetlabL : "GelPurify {{step.input[0]}}",
			template_wetlabR : "({{step.input[1] || 'L'}}, {{step.output}})"
		}
	});