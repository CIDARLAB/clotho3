angular.module('clotho.construction')
/**
 * @ngdoc value
 * @name ConstructionReactions
 * @description map of reactions to Clotho Function IDs and other info about them. Add to this list as new ones are supported
 */
	.value('ConstructionReactions', {
		"pcr": {
			reaction: "",
			clientReaction: "PCR.predict",
			readable: "PCR",
			template : "views/_construction/step/pcr.html"
		},
		"ligate": {
			reaction: "",
			clientReaction: "PCR.ligate",
			readable: "Ligation",
			template : "views/_construction/step/ligate.html"
		},
		"digest": {
			reaction: "",
			clientReaction: "Digest.digest",
			readable: "Restriction Digest",
			template : "views/_construction/step/digest.html"
		},
		"gelpurify": {
			reaction: "",
			clientReaction: "Digest.gelPurify",
			readable: "Gel Purify",
			template : "views/_construction/step/gelpurify.html"
		}
	});