angular.module('clotho.construction')
/**
 * @ngdoc value
 * @name ConstructionReactions
 * @description map of reactions to Clotho Function IDs and other info about them. Add to this list as new ones are supported
 */
	.value('ConstructionReactions', {
		"pcr": {
			reactionId: "clotho.functions.dna.pcrNucseq",
			clientReaction: "PCR.predict",
			readable: "PCR",
			template : "views/_construction/step/pcr.html"
		},
		"ligate": {
			reactionId: "clotho.functions.dna.ligate",
			clientReaction: "PCR.ligate",
			readable: "Ligation",
			template : "views/_construction/step/ligate.html"
		},
		"digest": {
			reactionId: "clotho.functions.dna.digest",
			clientReaction: "Digest.digest",
			readable: "Restriction Digest",
			template : "views/_construction/step/digest.html"
		},
		"gelPurify": {
			reactionId: "clotho.functions.dna.gelPurify",
			clientReaction: "Digest.gelPurify",
			readable: "Gel Purify",
			template : "views/_construction/step/gelpurify.html"
		}
	});