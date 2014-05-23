angular.module('clotho.construction')
/**
 * @ngdoc Value
 * @name constructionReactions
 * @description map of reactions to Clotho Function IDs and other info about them
 */
	.value('ConstructionReactions', {
		"pcr": {
			reaction: "PCR.predict",
			readable: "PCR",
			template : "views/_construction/step/pcr.html"
		},
		"ligate": {
			reaction: "PCR.ligate",
			readable: "Ligation",
			template : "views/_construction/step/ligate.html"
		},
		"digest": {
			reaction: "Digest.digest",
			readable: "Restriction Digest",
			template : "views/_construction/step/digest.html"
		},
		"gelpurify": {
			reaction: "Digest.gelPurify",
			readable: "Gel Purify",
			template : "views/_construction/step/gelpurify.html"
		}
	});